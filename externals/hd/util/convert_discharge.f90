! convert_discharge.f90 - Convert HD discharge to preselected ocean grids 
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

    PROGRAM convert_discharge
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
!     ***** Version 2.0 - June 2020
!           Programmierung und Entwicklung: Stefan Hagemann 
!           Generalization and modularisation of hdtonemo.f90
!           Input file:   rivmouth_source.nc
!                         coast_oceangrid.nc     1: Mouth, 0: others  (opposite to previous def.)
!           Namelist:     convert_inflow_ctl.nml
!           Output file:  hdcouple_<source model>_to_<ocean model>_imode<mode>.nc
!                         e.g. hdcouple_vs5_to_nemo_imode2.nc
!
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
    use mo_interpol,     ONLY: mapping, generate_point_mapping, write_mapping_data, &
                               define_separation, island_mapping
    use mo_convert,      ONLY: convert_inflow
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
    TYPE(model)          :: mask_target ! Mask with inflow points and target coordinates
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
    REAL, DIMENSION(:,:), ALLOCATABLE :: fdat
!
!   *** NETCDF variables
    INTEGER :: ierr, ncid, varid
    CHARACTER (len=20) :: clon, clat, cdimlon, cdimlat
!  
    CHARACTER (len=128) :: order, dnout
    CHARACTER (LEN=128) :: ctitle  ! Description of converted data
    INTEGER :: ipos
    CHARACTER (len=2) :: cmode     ! Command line input parameter for program mode  
    INTEGER :: imode               !  --> Method of using potential mouth masks
                   ! 1    Use existing mask with potential mouth points on ocean grid,e.g. coast_ocean_NEMO.nc
                   ! 2    Generate mask of coastal ocean points from ocean sea mask (Def.)
                   ! 3    Combine methods 1 and 2
                   ! 4    As 3, but with both masks prescribed
    CHARACTER (len=2) :: cin     ! Command line input parameter for source mouth mask 
    INTEGER :: isrc    = 2        ! Source data ID: 
                                 !   1    HD model Vs. 4
                                 !   2    HD model Vs. 5
                                 !   3    HD model Vs. 1.11
                                 !   4    mHm
                                 !   5    Utes standard input
                                 !   6,7    MPIOM
    CHARACTER (len=2) :: cocean  ! Command line input parameter for ID iocean   
    INTEGER :: iocean = 1        ! Ocean model setup ID: 
                                 !   1    NEMO Vs. 3.3/3.6 - North & Baltic Seas
                                 !   2    ECOSMO Vs. 3 - North & Baltic Seas
                                 !   3    SCHISM
                                 !   4    ECOSMO Vs. II - North & Baltic Seas
                                 !   5    ICON
                                 !   6    Nils
                                 !   7    ICON-Coast
                                 !   8    NEMO Vs. 4.0- North & Baltic Seas
                                 !   9    TRIM             14 MOM
!
!   *** Parameters to steer the conversion of discharge/bgc inflow: convert_inflow_ctl
    CHARACTER (len=192) :: dn_inflow = "${EXP}_meanflow_YYYY.nc"
    CHARACTER (len=192) :: dn_outflow = "${CFLOW}_${CDIS}_on_${OM}${OVS}_YYYY.nc"      
    INTEGER :: iconv = 0    ! Method of conversion: 0 = no, 1=Daily data (def.),
                            ! 2=monthly clim., 3=5 bgc sequence
                            ! 4=General sequence of bgc variables and timesteps
    INTEGER :: isep_in      ! Input data in separate files (No/Yes=0/1)
    INTEGER :: isep_out = 1 ! Output data in separate files (No/Yes=0/1)
    INTEGER :: ybeg         ! Start year of data
    INTEGER :: yend         ! End year of data
    INTEGER :: lu_nml = 10
    NAMELIST /convert_inflow_ctl/ dn_inflow, dn_outflow, iconv, isep_in, isep_out, ybeg, yend
!
!   *** Read command line input
    CALL GETARG(1, cmode)    ! Mode (former inemou)
    CALL GETARG(2, cin)      ! Source , e.g. HD or mHm
    CALL GETARG(3, cocean)   ! Target, i.e. ocean model
    CALL check_names(cmode, cin, cocean, imode, isrc, iocean, mask_src%name, mask_target%name)
!
!   *** Input files
    dn_src="rivmouth_source.nc"
    dn_target="coast_oceangrid.nc"
!
!   ******** Read grid info on source grid
    cdimlon = 'lon' ; cdimlat = 'lat'
    clon = 'lon' ; clat = 'lat'
    idim = 1              ! 1 coordinate dimension(s) (HD Parafile has 1)
    IF (isrc.EQ.4) THEN   ! mHm
      idim = 2            ! 2 coordinate dimensions
    ELSE IF (isrc.EQ.5) THEN    ! Utes data
      cdimlon = 'lon_dim' ; cdimlat = '-'
    ELSE IF (isrc.EQ.6) THEN    ! MPIOM BGC data
      cdimlon = 'r' ; cdimlat = 'c'
      clon = 'lon_2' ; clat = 'lat_2'
    ELSE IF (isrc.EQ.7) THEN    ! MPIOM data
      cdimlon = 'x' ; cdimlat = 'y'
    ENDIF
    WRITE(*,*) 'Source Grid No. ', isrc
    CALL read_grid_info(dn_src, gr_src, xmiss, clon, clat, cdimlon,cdimlat, idim)
    mask_src%nlon = gr_src%nlon
    mask_src%nlat = gr_src%nlat
!   *** Set specific values if, e.g., no grid is used.
    IF (isrc.EQ.5) THEN    ! Utes data
      gr_src%resolution = 0.5
      gr_src%origin_lon = -180.  ; gr_src%origin_lon = 90.
      gr_src%kshift_s = 0  ; gr_src%kshift_n = 0
      gr_src%kshift_w = 0  ; gr_src%kshift_e = 0
    ENDIF
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
    ierr = nf90_get_var(ncid,varid,mask_src%value)
    WRITE(*,*) 'Mouth array on source grid: ', &
                MINVAL(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps), ' - ', MAXVAL(mask_src%value)
    WRITE(*,*) '         SUM: ', SUM(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps)
!	  
!   ******** Read grid info on target grid, Mouth points and coordinates
!
!   *** 2D grid or ICON grid
    cdimlon = 'lon' ; cdimlat = 'lat'
    clon = 'lon' ; clat = 'lat'
    idim = 2                ! 2 dimensions
    IF (iocean.EQ.5 .OR. iocean.EQ.7) THEN    ! ICON grid --> NX=ncells, NY=1
      idim = 1              ! 1 dimensions
      cdimlon = 'ncells' ; cdimlat = '-'
      clon = 'clon' ; clat = 'clat'
    ELSE IF (iocean.EQ.3) THEN  ! SCHISM
      idim = 1              ! 1 dimension
      cdimlon = 'ncells' ; cdimlat = '-'
    ELSE IF (iocean.EQ.6) THEN
      idim = 1              ! 1 dimension
    ELSE IF (iocean.LE.4) THEN   ! NEMO 3.6/ECOSMO
      cdimlon = 'x' ; cdimlat = 'y'
    ELSE IF (iocean.LE.9) THEN   ! TRIM
      cdimlon = 'nx' ; cdimlat = 'ny'
    ENDIF
    WRITE(*,*) " =================================================================== "
    WRITE(*,*) 'Target Grid No. ', iocean, TRIM(dn_target)
    CALL read_grid_info(dn_target, gr_target, xmiss, clon, clat, cdimlon,cdimlat, idim)
    mask_target%nlon = gr_target%nlon
    mask_target%nlat = gr_target%nlat
    WRITE(*,*) "Target ocean model:  NX =" , gr_target%nlon, "  NY = ", gr_target%nlat
!
!   *** Feld Dimensionierung
    ALLOCATE(mask_target%xlon(mask_target%nlon,mask_target%nlat))
    ALLOCATE(mask_target%xlat(mask_target%nlon,mask_target%nlat))
    ALLOCATE(mask_target%value(mask_target%nlon,mask_target%nlat))
    ALLOCATE(fdat(mask_target%nlon,mask_target%nlat))

    IF (imode.EQ.3 .OR. imode.EQ.4) ALLOCATE(mask_target_preset(gr_target%nlon,gr_target%nlat))
    ALLOCATE(mask_on_ocean(gr_target%nlon,gr_target%nlat))
!
!   ******** READ mouth mask on source grid, e.g. HD: mask_src%value
    CALL read_coordinates(dn_target, mask_target, clon, clat)
    IF (iocean.EQ.5) THEN    ! ICON grid --> NX=ncells, NY=1
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
    ierr = nf90_open(dn_target, NF90_NOWRITE, ncid)
    ierr = nf90_inq_varid(ncid,'mask_coast_ocean',varid)
    ierr = nf90_get_var(ncid,varid,fdat)
    mask_target%value = NINT(fdat)
    WRITE(*,*) 'Mouth array on ocean grid: ', &
                MINVAL(mask_target%value), ' - ', MAXVAL(mask_target%value)
    WRITE(*,*) '         SUM: ', SUM(mask_target%value)

    IF (imode.EQ.3 .OR. imode.EQ.4) THEN
      map_char%nmask=2
      ierr = nf90_inq_varid(ncid,'nemo_mask_preset',varid)
      ierr = nf90_get_var(ncid,varid,fdat)
      mask_target_preset = NINT(fdat)
      WRITE(*,*) 'No. of preset primary mouths: ', SUM(mask_target_preset)
    ENDIF
    DEALLOCATE(fdat)
!
!   *** Define mapping characeristics search as maximum search radius
    CALL set_mapping_char(isrc, iocean, gr_src%resolution, map_char)
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
!
!   ********* Correction for Specific grid boxes
    CALL CORROCEAN(isrc, iocean, gr_src%nlon, gr_src%nlat, mask_src_mapped, &
         ix_target, iy_target, gr_target%nlon, gr_target%nlat, mask_on_ocean)

    WRITE(*,*) 'Maximum sources an ocean model gridbox is sharing: ', MAXVAL(mask_src_mapped)
!
!   ******** Schreiben der Index Arrays with ocean model River-Mouth targets
!
!   *** WRITE NETCDF output
    dnout = 'hdcouple_' // TRIM(mask_src%name) // '_to_' // TRIM(mask_target%name) // &
            '_imode' // TRIM(cmode) // '.nc'
    CALL write_mapping_data(dnout, mask_src, mask_target, imode, &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)
!
!   *** Convert discharge or bgc inflows
!
!   *** Read Namelist convert_inflow_ctl.nml
    INQUIRE(FILE='convert_inflow_ctl.nml', EXIST=logque)
    IF (logque) THEN
      OPEN (UNIT=lu_nml, FILE='convert_inflow_ctl.nml', STATUS='OLD')
      READ (lu_nml, NML=convert_inflow_ctl)
      CLOSE (lu_nml)
    ENDIF
    IF (iconv.gt.0) THEN
      ipos = INDEX(dn_inflow, '/', BACK=.TRUE.)
      ctitle = 'Conversion of ' // TRIM(dn_inflow(ipos+1:)) // &
               ' from ' // TRIM(mask_src%name) // ' to ' // TRIM(mask_target%name)
      CALL convert_inflow(dn_inflow, dn_outflow, ctitle, iconv, ybeg, yend, isep_in, isep_out, &
             mask_src, mask_target, xmiss, &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)      
    ENDIF
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
!
!   *********************************************************************************** 
    SUBROUTINE check_names(cmode, cin, cocean, imode, isrc, iocean, csrc, ctarget)
!   *********************************************************************************** 

      CHARACTER (LEN=*), INTENT(in)  :: cmode         
      CHARACTER (LEN=*), INTENT(in)  :: cin         
      CHARACTER (LEN=*), INTENT(in)  :: cocean         
      INTEGER, INTENT(out)           :: imode     ! Mode ID, former inemou
      INTEGER, INTENT(out)           :: isrc       ! Source data ID
      INTEGER, INTENT(out)           :: iocean    ! Ocean model ID
      CHARACTER (LEN=*), INTENT(out)  :: csrc     ! Source name         
      CHARACTER (LEN=*), INTENT(out)  :: ctarget  ! Target name         
!
      INTEGER, PARAMETER :: iocmax = 15
      INTEGER, PARAMETER :: isrcmax = 7
      CHARACTER (len=20), DIMENSION(isrcmax) :: cname_src  ! Source model names
      CHARACTER (len=20), DIMENSION(iocmax) :: cname_oc   ! Ocean model names 
!
      DATA cname_src / 'hd_vs4', 'hd_vs5', 'hd_vs1_11', 'mhm_vs2',  'ute', 'mpiom-bgc', 'mpiom' / 
      DATA cname_oc / 'nemo', 'ecosmo3', 'SCHISM', 'ecosmo2', 'iconomip', &
                      'nsea', 'icon-coast', 'nemo_vs4', 'trim', 'nemo-med7km', &
                      'hd_vs4', 'hd_vs5', 'hd_vs1_11', 'mom', 'nemo_nss'  / 
!
      READ(cmode, '(I2)') imode
      IF (imode.LT.1 .OR. imode.GT.4) THEN
        WRITE(*,*) ' imode out of range [1,4]: ', imode
        STOP 'ERROR in check_names --> TERMINATED!'
      ENDIF
!
      READ(cin, '(I2)') isrc
      IF (isrc.LT.1 .OR. isrc.GT.isrcmax) THEN
        WRITE(*,'(A,I2,A,I2)') ' isrc out of range [1,',isrcmax, ']: ', isrc
        STOP 'ERROR in check_names --> TERMINATED!'
      ENDIF
      csrc = cname_src(isrc)

      READ(cocean, '(I2)') iocean
      IF (iocean.LT.1 .OR. iocean.GT.iocmax) THEN
        WRITE(*,'(A,I2,A,I2)') ' iocean out of range [1,',iocmax, ']: ', iocean
        STOP 'ERROR in check_names --> TERMINATED!'
      ENDIF
      ctarget = cname_oc(iocean)

    END SUBROUTINE check_names

!   *********************************************************************************** 
      SUBROUTINE CORROCEAN(isrc, iocean, NL, NB, mask_src_mapped, ix_target, iy_target, NX, NY, mask_on_ocean)
!   *********************************************************************************** 
!
!     ***  Correction of specific grid boxes depending on the ocean model grid.
!     ***  1) Removal of inflow points, e.g. near the ocean model domain boundary 
!     ***  2) Choosing a dedicated inflow point, e.g. for the Elbe
! 
      INTEGER, INTENT(in) :: isrc      ! Source model setup ID
      INTEGER, INTENT(in) :: iocean    ! Ocean model setup ID
      INTEGER, INTENT(in) :: NL, NB
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: mask_src_mapped  ! Mask with HD mouths with a valid 
                                                            ! nearest mouth on ocean grid
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: ix_target    ! x-Indices of nearest ocean grid mouth points
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: iy_target    ! y-Indices of nearest ocean grid mouth points

      INTEGER, INTENT(in) :: NX, NY
      INTEGER, DIMENSION(NX,NY), INTENT(inout) :: mask_on_ocean  ! Mask with associated ocean mouth points

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

      SELECT CASE (iocean)
        CASE (1)
          IF (NL.EQ.960 .AND. NB.EQ.540) THEN 
            NCORR=1
            ALLOCATE(IX(NCORR)) ; ALLOCATE(IY(NCORR))
            ALLOCATE(INEWX(NCORR)) ; ALLOCATE(INEWY(NCORR))
            IX(1) = 495  ; IY(1) = 145 ; INEWX(1) = 901 ; INEWY(1) = 598  ! Newa -> NEMO sets Eastern boundary to NAN
          ENDIF
        CASE (2)

!!          IF (NL.EQ.832 .AND. NB.EQ.592) THEN        ! mRm source
          IF (isrc.EQ.4) THEN        ! mRm source
            NCORR=2
            ALLOCATE(IX(NCORR)) ; ALLOCATE(IY(NCORR))
            ALLOCATE(INEWX(NCORR)) ; ALLOCATE(INEWY(NCORR))
            IX(1) = 262  ; IY(1) = 324 ; INEWX(1) = 248 ; INEWY(1) = 709  ! Maass -> Separation from Rhine
            IX(2) = 250  ; IY(2) = 325 ; INEWX(2) = 248 ; INEWY(2) = 701  ! Put Rhine from 709 analogous to HD 
            CALL define_separation(IX(1), IY(1), 250, 325)                ! Separation of Maass from Rhine
            mask_src%value(IX(1), IY(1)) = 1

          ELSE IF (NL.EQ.960 .AND. NB.EQ.540) THEN   ! HD Euro 5 Min. source
            NREM=2
            ALLOCATE(IREMX(NREM)) ; ALLOCATE(IREMY(NREM))
            IREMX(1) = 281  ; IREMY(1) = 334
            IREMX(2) = 294  ; IREMY(2) = 341
          ENDIF
        CASE (4)
!!          IF (NL.EQ.832 .AND. NB.EQ.592) THEN        ! mRm source
          IF (isrc.EQ.4) THEN        ! mRm source
            NCORR=6
            ALLOCATE(IX(NCORR)) ; ALLOCATE(IY(NCORR))
            ALLOCATE(INEWX(NCORR)) ; ALLOCATE(INEWY(NCORR))
!
!           *** Removal mRm
            NREM=22
            ALLOCATE(IREMX(NREM)) ; ALLOCATE(IREMY(NREM))
            IREMX(1) = 11  ; IREMY(1) = 69   ! Orkney
            IREMX(2) = 11  ; IREMY(2) = 70   ! Orkney
            IREMX(3) = 11  ; IREMY(3) = 71   ! Orkney
            IREMX(4) = 10  ; IREMY(4) = 72   ! Orkney
            IREMX(5) = 9   ; IREMY(5) = 73   ! Orkney
            IREMX(6) = 9   ; IREMY(6) = 74   ! Orkney
            IREMX(7) = 14  ; IREMY(7) = 66   ! Shetland
            IREMX(8) = 16  ; IREMY(8) = 66   ! Shetland
            IREMX(9) = 17  ; IREMY(9) = 66   ! Shetland
            IREMX(10) = 18 ; IREMY(10) = 66  ! Shetland
            IREMX(11) = 19 ; IREMY(11) = 66  ! Shetland
            IREMX(12) = 56 ; IREMY(12) = 66  ! Norwegian coast
            IREMX(13) = 57 ; IREMY(13) = 66  ! Norwegian coast
            IREMX(14) = 57 ; IREMY(14) = 67  ! Norwegian coast
            IREMX(15) = 57 ; IREMY(15) = 68  ! Norwegian coast
            IREMX(16) = 57 ; IREMY(16) = 69  ! Norwegian coast
            IREMX(17) = 16 ; IREMY(17) = 156 ! UK coast (too close to boundary)
            IREMX(18) = 17 ; IREMY(18) = 156 ! UK coast (too close to boundary)
            IREMX(19) = 16 ; IREMY(19) = 164 ! French coast (too close to boundary)
            IREMX(20) = 17 ; IREMY(20) = 164 ! French coast (too close to boundary)
            IREMX(21) = 19 ; IREMY(21) = 166 ! French coast (too much runoff)
            IREMX(22) = 20 ; IREMY(22) = 167 ! French coast (too much runoff)
!           *** Correction
            IX(1) = 334  ; IY(1) = 296  ;  INEWX(1) = 78  ; INEWY(1) = 120    ! Elbe
            IX(2) = 270  ; IY(2) = 311  ;  INEWX(2) = 57  ; INEWY(2) = 128    ! Ijssel
            IX(3) = 525  ; IY(3) = 101  ;  INEWX(3) = 158 ; INEWY(3) = 7      ! Lule
            IX(4) = 514  ; IY(4) = 274  ;  INEWX(4) = 145 ; INEWY(4) = 115    ! Pregel
            IX(5) = 250  ; IY(5) = 325  ;  INEWX(5) = 50 ; INEWY(5) = 141     ! Rhein -> same target as for HD
            IX(6) = 262  ; IY(6) = 324  ;  INEWX(6) = 50 ; INEWY(6) = 142     ! Maass -> Separation from Rhine
            CALL define_separation(IX(6), IY(6), 250, 325)                    ! Separation of Maass from Rhine
            mask_src%value(IX(6), IY(6)) = 1
          ELSE IF (isrc.EQ.1 .OR. isrc.EQ.2) THEN   ! HD Vs. 4 & 5 - Euro 5 Min. source
            NCORR=9
            ALLOCATE(IX(NCORR)) ; ALLOCATE(IY(NCORR))
            ALLOCATE(INEWX(NCORR)) ; ALLOCATE(INEWY(NCORR))
!
!           *** Removal HD - HD 4 and 5
            NREM=8
            ALLOCATE(IREMX(NREM)) ; ALLOCATE(IREMY(NREM))
            IREMX(1) = 18  ; IREMY(1) = 66   ! Shetland
            IREMX(2) = 57  ; IREMY(2) = 66   ! Norwegian coast
            IREMX(3) = 9   ; IREMY(3) = 73   ! Orkney
            IREMX(4) = 16  ; IREMY(4) = 156  ! UK coast
            IREMX(5) = 16  ; IREMY(5) = 164  ! French coast (too close to boundary)
            IREMX(6) = 18  ; IREMY(6) = 164  ! French coast (too close to boundary)
            IREMX(7) = 19  ; IREMY(7) = 166  ! French coast (too much runoff)
            IREMX(8) = 20  ; IREMY(8) = 167  ! French coast (too much runoff)
!
!           *** Correction
            IF (isrc.EQ.1) THEN  
              IX(1) = 251  ; IY(1) = 222  ;  INEWX(1) = 78  ; INEWY(1) = 120    ! Elbe
            ELSE IF (isrc.EQ.2) THEN  
              IX(1) = 246  ; IY(1) = 220  ;  INEWX(1) = 78  ; INEWY(1) = 120    ! Elbe
            ENDIF
            IX(2) = 200  ; IY(2) = 233  ;  INEWX(2) = 57  ; INEWY(2) = 128    ! Ijssel (before 54, 134)
            IX(3) = 393  ; IY(3) =  75  ;  INEWX(3) = 158 ; INEWY(3) = 7      ! Lule   (before 161, 4)
            IX(4) =  81  ; IY(4) = 175  ;  INEWX(4) = 5   ; INEWY(4) = 83     ! Ness
            IX(5) =  79  ; IY(5) = 175  ;  INEWX(5) = 5   ; INEWY(5) = 83     ! Beauly
            IX(6) =  79  ; IY(6) = 174  ;  INEWX(6) = 5   ; INEWY(6) = 83     ! Conon
            IX(7) =  80  ; IY(7) = 170  ;  INEWX(7) = 5   ; INEWY(7) = 82     ! Kyle of Sutherland
            IX(8) = 114  ; IY(8) = 255  ;  INEWX(8) = 18  ; INEWY(8) = 156    ! UK Southcoast (before 17,156)
            IX(9) = 114  ; IY(9) = 256  ;  INEWX(9) = 18  ; INEWY(9) = 156    ! UK Southcoast (before 17,156)
          ENDIF
        CASE DEFAULT
           RETURN
      END SELECT
!
!     *** Removal of ocean mouths
      IF (NREM.GE.1) THEN
      DO I = 1, NREM
        NDUM = mask_on_ocean(IREMX(I), IREMY(I))
        DO JB=1, NB
        DO JL=1, NL
        IF (ix_target(JL,JB).EQ.IREMX(I) .AND. iy_target(JL,JB).EQ.IREMY(I)) THEN
          ix_target(JL,JB) = 0
          iy_target(JL,JB) = 0
          mask_src_mapped(JL,JB) = 0
          mask_on_ocean(IREMX(I), IREMY(I)) = mask_on_ocean(IREMX(I), IREMY(I)) - 1
        ENDIF
        ENDDO
        ENDDO
        WRITE(*,*) I, '. point reduced from ', NDUM, ' to ', mask_on_ocean(IREMX(I), IREMY(I))
      ENDDO
      DEALLOCATE(IREMX)
      DEALLOCATE(IREMY)
      ENDIF
!
!     *** Correction of target points
      IF (NCORR.GE.1) THEN
        DO I = 1, NCORR
          JL = IX(I) ; JB = IY(I)
          mask_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) = mask_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) - 1
          ix_target(JL,JB) = INEWX(I)
          iy_target(JL,JB) = INEWY(I)
          mask_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) = mask_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) + 1
          mask_src_mapped(JL,JB) = mask_on_ocean(ix_target(JL,JB), iy_target(JL,JB))
          WRITE(*,*) I, '. point at ', JL,JB, ' corrected towards ', ix_target(JL,JB), iy_target(JL,JB)
        ENDDO
        DEALLOCATE(IX) ; DEALLOCATE(IY)
        DEALLOCATE(INEWX) ; DEALLOCATE(INEWY)
      ENDIF

      END SUBROUTINE CORROCEAN

!   *********************************************************************************** 
    SUBROUTINE set_mapping_char(isrc, iocean, res_src, map_char)
!   *********************************************************************************** 

      INTEGER, INTENT(in) :: isrc  ! Source data model ID: 
                                   !   1    HD model Vs. 4
                                   !   2    HD model Vs. 5
                                   !   3    HD model Vs. 1.10
                                   !   4    mHm
                                   !   5    Utes standard input
                                   !   6    MPIOM-BGC
                                   !   7    MPIOM
      INTEGER, INTENT(in) :: iocean  ! Target Ocean model setup ID: 
                                   !   1    NEMO Vs. 3.3/3.6 - North & Baltic Seas
                                   !   2    ECOSMO Vs. 3 - North & Baltic Seas
                                   !   3    SCHISM
                                   !   4    ECOSMO Vs. II - North & Baltic Seas
                                   !   5    ICON
                                   !   6    Nils Nordseemodell
                                   !   7    ICON-Coast
                                   !   8    NEMO Vs. 4.0 - North & Baltic Seas
                                   !   9    TRIM - North & Baltic Seas
                                   !  10    NEMO-med7km - Mediterranean Sea
                                   !  11    HD model Vs. 4
                                   !  12    HD model Vs. 5
                                   !  13    HD model Vs. 1.10
                                   !  14    MOM - Baltic Sea
                                   !  15    NEMO-NSS - North Sea
      DOUBLE PRECISION, INTENT(in) :: res_src    ! Source grid resolution, e.g. HD model
      TYPE(mapping), INTENT(out) :: map_char     ! Mapping charactistics

      map_char%deg_max_bound=res_src       ! 1 Gridbox outside ocean domain
      IF (isrc.EQ.3) THEN         ! 0.5 degree
        map_char%dist_max_prime=25000.     ! 25 km for primary mask
        map_char%dist_max=100000.          ! 100 km for secondary mask (default if 1 mask)
        IF (iocean.EQ.5) map_char%dist_max=400000.    ! ICON ocean 60km-coast --> 200 km 
                                           ! For all bgc inflows, 1000 km necessary
        IF (iocean.EQ.7) THEN 
          map_char%dist_max_prime=60000.   ! 
          map_char%dist_max=600000.        ! 
        ENDIF
      ELSE IF (isrc.LE.2) THEN    ! 5 Min.
        map_char%dist_max_prime=4000.     ! 4 km for primary mask
        map_char%dist_max=16000.    ! 16 km for secondary mask (default for only 1 mask)
        IF (iocean.EQ.1 .OR. iocean.EQ.8) map_char%dist_max=200000.  ! NEMO ocean coast too smooth --> 200 km 
        IF (iocean.EQ.2) map_char%dist_max=100000.
        IF (iocean.EQ.3) map_char%dist_max=40000.   ! Necessary in Northern Russia 
        IF (iocean.EQ.4) map_char%dist_max=100000.  ! ECOSMO-10 km ocean coast very smooth --> 200 km 
        IF (iocean.EQ.6) map_char%dist_max=57000.   ! N Sea model: Rhine are in, Baltc Sea out
        IF (iocean.EQ.9) map_char%dist_max=200000.  ! TRIM ocean coast is smooth --> 200 km 
        IF (iocean.EQ.10) map_char%dist_max=80000.  ! Necessary in Greece 
        IF (iocean.EQ.14) map_char%dist_max=50000.  ! 
        IF (iocean.EQ.15) map_char%dist_max=67000.  ! 
      ELSE IF (isrc.EQ.4) THEN   ! mHm
        map_char%dist_max_prime=4000.               ! 4 km for primary mask
        map_char%dist_max=16000.        ! 16 km for secondary mask (default for only 1 mask)
        map_char%deg_max_bound=0.5
        IF (iocean.EQ.1 .OR. iocean.EQ.8) map_char%dist_max=200000.  ! NEMO ocean coast too smooth --> 200 km 
        IF (iocean.EQ.2) map_char%dist_max=100000.  ! UFZ data contain sinks that should be 
                                                    ! connected to the coast 
        IF (iocean.EQ.4) map_char%dist_max=100000.  ! ECOSMO-10 km ocean coast very smooth --> 200 km 
      ELSE IF (isrc.EQ.6) THEN   ! MPIOM
        map_char%dist_max_prime=60000.     ! 60 km for primary mask
        map_char%dist_max=200000.          ! 200 km for secondary mask (default if 1 mask)
      ELSE IF (isrc.EQ.7) THEN   ! MPIOM
        map_char%dist_max_prime=60000.     ! 60 km for primary mask
        map_char%dist_max=200000.          ! 200 km for secondary mask (default if 1 mask)
      ENDIF
      WRITE(*,*) "Ocean model ", icocean, " Inflow resolution: ", res_src   !  , ' =?', 0.5/6. 
      WRITE(*,*) "Maximum distances: Primary = ", map_char%dist_max_prime 
      WRITE(*,*) "                 Secondary = ", map_char%dist_max 
      WRITE(*,*) "                  Boundary = ", map_char%deg_max_bound 
!
    END SUBROUTINE set_mapping_char
!	  
    END PROGRAM convert_discharge
