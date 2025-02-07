! hdtocoean.f90 - Links the HD river mouths to coastal inflow points on a given ocean grid.
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM HDTOOCEAN
!
!     ***** Programm that takes the HD river mouth points that have been created from the
!           River-Direction-File (HD para Array FDIR) beforehand, either on the HD grid
!           or already transformed to the NEMO grid. Then it finds the 
!           the nearest river mouth points on the NEMO grid. 
!           The latter is read in from file into mask_target%value. The program 
!           creates an HD array with the NEMO lon and lat indices that
!           allocate these nearest mouth points to each HD mouth point.
!           In addition, discharge/ bgc inflow may be converted using the
!           generated mapping data. This is steered via the namelist convert_inflow_ctl.nml
!
!     ***** Output
!     File:     hd_to_ocean_mouth.nc
!     Variable: FMOU_HD_TO_NEMO   HD Mouth point with an associated ocean model, e.g. NEMO, mouth point
!                                 on the source grid (HD or ocean model grid)
!               INDEXX            Variable ix_target, Ocean model longitude index 
!                                 of associated mouth for the HD mouth point
!               INDEXY            Variable iy_target, Ocean model latitude index 
!                                 of associated mouth for the HD mouth point
!       nemo_hdmouth.nc  --> hdmouth_on_oceangrid.nc
!           Variable: FMOU_HD_ON_NEMO   Mask of associated NEMO mouth points on NEMO grid
!
!     ***** Version 2.1 - January 2022
!           Programmierung und Entwicklung: Stefan Hagemann 
!           Separation and generalization of routines from convert_discharge.f90
!           Input file:   rivmouth_source.nc
!                         coast_oceangrid.nc     1: Mouth, 0: others
!           Namelist:     hdtoocean_ctl.nml
!           Output file:  hdcouple_<Experiment ID/TAG>.nc
!                         e.g. hdcouple_vs5_to_nemo.nc
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
!     ***        Anmerkung: Richtung 5 = Discharge-Trap = Sink on land
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
!     **** Variablenliste:
!     ***  
!     ***  dn_src = Name des RDF
!     ***   ique = Kommentarvariable ( 0 = Kein Kommentar )
!     *** FMOUTH = Muendungsmaske HD
!     ***   mask_target%value = Muendungsmaske NEMO
!     ***   FDIR = HD Direction array (not used)
!     *** ix_target = x-Indices of nearest NEMO mouth points
!     *** iy_target = y-Indices of nearest NEMO mouth points
!     ***  
!     *** DISMAX = Distance threshold
!     *** DEGMAX = Maximum Distance in degree from outside the NEMO region
!     ***
    use netcdf
    use mo_grid,         ONLY: domain, model, read_grid_info, read_coordinates
    use mo_interpol,     ONLY: mapping, generate_point_mapping, island_mapping, write_mapping_data
!
    DOUBLE PRECISION, PARAMETER :: xmiss=-9999. 
    DOUBLE PRECISION, PARAMETER :: PI=3.14159265
    DOUBLE PRECISION, PARAMETER :: zeps = 1.E-6 
 
    INTEGER, PARAMETER  :: ique =0     ! Print out some debug statements for ique=1
!
!   *** Define source and target grid
!
!   *** Grid info
    CHARACTER (len=128) :: dn_src      ! File with inflow mask on source grid
    CHARACTER (len=128) :: dn_target   ! File with potential inflow mask on target grid
    TYPE(domain)        :: gr_src      ! Source grid
    TYPE(domain)        :: gr_target   ! Target grid
    TYPE(model)         :: mask_src    ! Mask with inflow points and source coordinates
    TYPE(model)         :: mask_target ! Mask with inflow points and target coordinates
    TYPE(mapping)       :: map_char    ! Mapping characteristics

    INTEGER, DIMENSION(2) :: ICOORD
!
!   *** Arrays on source grid, e.g. HD
!
!   Mask of source mouths with a valid nearest target mouth
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_src_mapped   
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ix_target
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: iy_target
!   Ocean model (e.g. NEMO) arrays
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_target_preset ! Note: Regular 0/1 mask
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_on_ocean
!
!   *** NETCDF variables
    INTEGER :: ierr, ncid, varid
    CHARACTER (len=20) :: clon, clat, cdimlon, cdimlat
    INTEGER :: ivar_type
!  
    REAL, DIMENSION(:,:), ALLOCATABLE :: fdum
    CHARACTER (len=128) :: order, dnout
    INTEGER, PARAMETER :: imode = 2 ! Use mask of coastal ocean points derived from ocean sea mask
!
!   *** Parameters to steer the conversion of discharge/bgc inflow: convert_inflow_ctl
    CHARACTER (len=80)  :: cexp = "hd_to_ocean"     ! Experiment ID
    INTEGER :: irad = 0            ! 0=no trafo, 1 = coordinate trafo from radian to degree
    INTEGER :: icon = 0            ! 0=no ICON grid, 1 = ICON grid, i.e. only 1 coordinate dimension
    DOUBLE PRECISION :: search_radius = 0.
    DOUBLE PRECISION :: search_boundary = 0.
    LOGICAL :: lexbound = .false.  ! Exclude boundary boxes for target grid             
    CHARACTER (len=128) :: name_src = 'HD'              ! Name of source grid/model
    CHARACTER (len=128) :: name_target = 'Ocean model'  ! Name of target grid/model
    NAMELIST /hdtoocean_ctl/ cexp, irad, icon, search_radius, search_boundary, lexbound, &
                             name_src, name_target  
    INTEGER :: lu_nml = 10
!
!   *** Input files
    dn_src="rivmouth_source.nc"
    dn_target="coast_oceangrid.nc"
!
!   *** Define mapping characeristics and read Namelist hdtoocean_ctl.nml
    INQUIRE(FILE='hdtoocean_ctl.nml', EXIST=logque)
    IF (logque) THEN
      OPEN (UNIT=lu_nml, FILE='hdtoocean_ctl.nml', STATUS='OLD')
      READ (lu_nml, NML=hdtoocean_ctl)
      CLOSE (lu_nml)
    ENDIF
    
    IF (search_radius.GT.zeps) THEN
      map_char%dist_max = search_radius * 1000.
    ELSE
      map_char%dist_max = 100000.
    ENDIF
    map_char%dist_max_prime = map_char%dist_max / 4.
    IF (search_boundary.GT.zeps) map_char%deg_max_bound = search_boundary
    WRITE(*,*) "  Maximum search radius = ", map_char%dist_max/1000., ' km' 
    WRITE(*,*) "Search outside boundary = ", map_char%deg_max_bound, ' degree' 
!
!   ******** Read grid info on HD source grid *****************************************
!
    cdimlon = 'lon' ; cdimlat = 'lat'
    clon = 'lon' ; clat = 'lat'
    idim = 1              ! 1 coordinate dimension(s) (HD Parafile has 1)
    WRITE(*,*) 'Source Grid: HD '
    CALL read_grid_info(dn_src, gr_src, xmiss, clon, clat, cdimlon,cdimlat, idim)
    mask_src%nlon = gr_src%nlon
    mask_src%nlat = gr_src%nlat
    mask_src%name = name_src
!
!   *** Feld Dimensionierung
    ALLOCATE(mask_src%value(gr_src%nlon,gr_src%nlat))
    ALLOCATE(mask_src%xlon(gr_src%nlon,gr_src%nlat))
    ALLOCATE(mask_src%xlat(gr_src%nlon,gr_src%nlat))
!
    ALLOCATE(mask_src_mapped(gr_src%nlon,gr_src%nlat))
    ALLOCATE(ix_target(gr_src%nlon,gr_src%nlat))
    ALLOCATE(iy_target(gr_src%nlon,gr_src%nlat))
!
!   ******** READ mouth mask on source grid, e.g. HD: mask_src%value
    CALL read_coordinates(dn_src, mask_src, clon, clat)
    WRITE(*,*) "Read Source mask --> Origin: ", mask_src%xlon(1,1), mask_src%xlat(1,1)
    WRITE(*,*) "                Lower right: ", &
         mask_src%xlon(gr_src%nlon,gr_src%nlat), mask_src%xlat(gr_src%nlon,gr_src%nlat)
    ierr = nf90_open(dn_src, NF90_NOWRITE, ncid)
    ierr = nf90_inq_varid(ncid,'FMOUTH',varid)
    ierr = nf90_inquire_variable(ncid, varid, xtype=ivar_type)
    WRITE(*,*) "              Variable Type: ", ivar_type
    IF (ivar_type.EQ.NF90_INT) THEN
      ierr = nf90_get_var(ncid,varid,mask_src_mapped)
      WRITE(*,*) 'ierr= ',ierr
      mask_src%value = mask_src_mapped
    ELSE IF (ivar_type.EQ.NF90_FLOAT) THEN
      ALLOCATE(fdum(gr_src%nlon,gr_src%nlat))
      ierr = nf90_get_var(ncid,varid,fdum)
      WRITE(*,*) 'ierr= ',ierr
      mask_src%value = int(fdum)
      DEALLOCATE(fdum)
    ELSE
       STOP 'Mask value is not an integer array -> STOP!'
    ENDIF
    WRITE(*,*) 'Mouth array on source grid: ', &
                MINVAL(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps), ' - ', MAXVAL(mask_src%value)
    WRITE(*,*) '         SUM: ', SUM(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps)
!	  
!   ******** Read grid info on target grid, Mouth points and coordinates *****************************
!
!   *** Default: 2D grid (ICON grid may work to)
    cdimlon = 'lon' ; cdimlat = 'lat'
    clon = 'lon' ; clat = 'lat'
    idim = 2                ! 2 dimensions
    IF (icon.EQ.1) THEN     ! ICON grid --> NX=ncells, NY=1
      idim = 1              ! 1 dimensions
      cdimlon = 'ncells' ; cdimlat = '-'
      clon = 'clon' ; clat = 'clat'
    ENDIF

    WRITE(*,*) " =================================================================== "
    WRITE(*,*) 'Read info of target ocean grid '
    CALL read_grid_info(dn_target, gr_target, xmiss, clon, clat, cdimlon,cdimlat, idim)
    mask_target%nlon = gr_target%nlon
    mask_target%nlat = gr_target%nlat
    mask_target%name = name_target
    WRITE(*,*) "Target ocean model:  NX =" , gr_target%nlon, "  NY = ", gr_target%nlat
!
!   *** Feld Dimensionierung
    ALLOCATE(mask_target%xlon(mask_target%nlon,mask_target%nlat))
    ALLOCATE(mask_target%xlat(mask_target%nlon,mask_target%nlat))
    ALLOCATE(mask_target%value(mask_target%nlon,mask_target%nlat))

    ALLOCATE(mask_on_ocean(gr_target%nlon,gr_target%nlat))
!
!   ******** READ mouth mask on source grid, e.g. HD: mask_src%value
    CALL read_coordinates(dn_target, mask_target, clon, clat)

    ! Coordinate trafo from radian to degree
    IF (irad.EQ.1) THEN    
      mask_target%xlon = mask_target%xlon / PI * 180.
      WHERE (mask_target%xlon.GT.180) mask_target%xlon = mask_target%xlon - 360.
      WHERE (mask_target%xlon.LT.-180) mask_target%xlon = mask_target%xlon + 360.
      mask_target%xlat = mask_target%xlat / PI * 180.
    ENDIF

    WRITE(*,*) "  Read target mask --> Origin: ", mask_target%xlon(1,1), mask_target%xlat(1,1)
    WRITE(*,*) "                  Lower right: ", &
         mask_target%xlon(mask_target%nlon,mask_target%nlat), &
         mask_target%xlat(mask_target%nlon,mask_target%nlat)
    WRITE(*,*) "Ocean Longitudes: ", MINVAL(mask_target%xlon), ' - ', MAXVAL(mask_target%xlon)
    WRITE(*,*) " Ocean Latitudes: ", MINVAL(mask_target%xlat), ' - ', MAXVAL(mask_target%xlat)
!     CALL read_variable(dn_src, cvar, mask_src%value)

    ! *** Read Target mouth point mask
    ! *** Note that 1 = Coastal Ocean Point, 0 = Other ocean points, 2 = land point

    WRITE(*,*) 'open: ', TRIM(dn_target)
    ierr = nf90_open(dn_target, NF90_NOWRITE, ncid)
    ierr = nf90_inq_varid(ncid,'mask_coast_ocean',varid)
    ierr = nf90_get_var(ncid,varid,mask_target%value)
    WRITE(*,*) '   Mouth array on ocean grid: ', &
                MINVAL(mask_target%value), ' - ', MAXVAL(mask_target%value)
    WRITE(*,*) ' No. of coastal ocean points:', COUNT(mask_target%value.EQ.1)
!
!   *** Exclude boundary boxes from target mask
    IF (lexbound) THEN
      WHERE(mask_target%value(1,:).EQ.1) mask_target%value(1,:) = 2 
      WHERE(mask_target%value(mask_target%nlon,:).EQ.1) mask_target%value(mask_target%nlon,:) = 2 
      WHERE(mask_target%value(:,1).EQ.1) mask_target%value(:,1) = 2 
      WHERE(mask_target%value(:,mask_target%nlat).EQ.1) mask_target%value(:,mask_target%nlat) = 2 
      WRITE(*,*) 'No. after excluding boundary: ', COUNT(mask_target%value.EQ.1)
    ENDIF
!
!   ******** Search for nearest target mouth point and generate Index arrays
    CALL generate_point_mapping(mask_src, mask_target, map_char,  &
              mask_target_preset, xmiss, ique,  &
              mask_src_mapped, ix_target, iy_target, mask_on_ocean)
!
!   *** If the ocean model resolution is coarser than source resolution, some islands on the source 
!   *** grid may not have coastal ocean points nearby. 
    CALL island_mapping(mask_src, mask_target, map_char, ique,  &
              mask_src_mapped, ix_target, iy_target, mask_on_ocean)

    WRITE(*,*) 'Source mouth points with unique ocean gridbox target: ', &
               SUM(mask_src_mapped, mask_src_mapped.EQ.1)
    WRITE(*,*) 'Source mouth points with shared ocean gridbox targets: ', &
               COUNT(mask_src_mapped.gt.1.5)

    WRITE(*,*) 'Maximum sources an ocean model gridbox is sharing: ', MAXVAL(mask_src_mapped)
!
!   ******** Schreiben der Index Arrays with ocean model River-Mouth targets
!
!   *** WRITE NETCDF output
    dnout = 'hdcouple_' // TRIM(cexp) // '.nc'
    CALL write_mapping_data(dnout, mask_src, mask_target, imode, &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)
!
!   *** The End
    DEALLOCATE(mask_src_mapped)
    DEALLOCATE(mask_src%value)
    DEALLOCATE(mask_src%xlon)
    DEALLOCATE(mask_src%xlat)
    DEALLOCATE(ix_target)
    DEALLOCATE(iy_target)
    DEALLOCATE(mask_target%value)
    IF (imode.EQ.3 .OR. imode.EQ.4) DEALLOCATE(mask_target_preset)
    DEALLOCATE(mask_target%xlon)
    DEALLOCATE(mask_target%xlat)
!
!     ******** Programmende
      CONTAINS


    END PROGRAM HDTOOCEAN
