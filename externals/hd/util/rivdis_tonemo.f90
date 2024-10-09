! rivdis_tonemo.f90 - Program that takes the remapping file and converts HD discharges into NEMO inflows
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM RIVDIS_TONEMO
!
!     ******** Programm that takes reads the HD river discharge at HD river mouth points 
!              and transfers it to the NEMO grid and the associated NEMO mouth points.
!              It reads the file  hd_nemo_mouth.nc and the arrays FMOU_HD_NEMO, INDEXX, INDEXY
!              NEMO grid info is read from the hdmouth_on_nemo.nc.
!
!     *** FMOU_HD_NEMO = Mask with HD mouth boxes that have an associated NEMO mouth point.
!     *** INDEXX = x-Indices of nearest NEMO mouth points
!     *** INDEXY = y-Indices of nearest NEMO mouth points
!
!              
!
!     ******** Version 1.0 - Dec. 2016
!              Programmierung und Entwicklung: Stefan Hagemann 
!
!     ******** Vs. 1.1 - Feb. 2017
!              Implement test output for a selected coordinate ILLOG, IBLOG
!              Elbe: Lon = 8.5, Lat=54.5   Cat=12
!                    HD:   378      72
!           HD_NEMO 3.3:   233     174   after nearest neighbor interpolation
!        NEMO 3.3-coast:   230     174
!           HD_NEMO 3.6:   517     426   after nearest neighbor interpolation (EHYPE)
!
      use netcdf
!
      PARAMETER(XMISS=-9999.) 
      INTEGER, PARAMETER  :: IQUE =0      ! Print out some debug statements for IQUE=1
      INTEGER, PARAMETER  :: ILLOG =378   ! Log Output for selected mouth point with lon index=ILLOG
      INTEGER, PARAMETER  :: IBLOG =72    ! Log Output for selected mouth point with lat index=IBLOG
      INTEGER NL, NB, NX, NY
!     HD arrays
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_MOU
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEXX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEXY
      REAL, DIMENSION(:,:), ALLOCATABLE :: FRIV
!     NEMO arrayS
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_NEMO
      REAL, DIMENSION(:,:), ALLOCATABLE :: XNEMO
      REAL, DIMENSION(:,:), ALLOCATABLE :: YNEMO
      REAL, DIMENSION(:,:), ALLOCATABLE :: FRIV_NEMO
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, dimid, varid, dimids(2), ncid2
!  
      CHARACTER (len=128) :: DNHDTONEMO, DNNEMO, ORDER, DNOUT, DNINP
!
      DNHDTONEMO="hd_to_nemo_mouth.nc"
      DNNEMO="hdmouth_on_nemo.nc"
      DNINP="discharge.nc"
!
!     ******** READ HD array MASK_MOU
      ierr = nf90_open(DNHDTONEMO, NF90_NOWRITE, ncid)

      ierr = nf90_inq_dimid(ncid,'lon',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NL)
      ierr = nf90_inq_dimid(ncid,'lat',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NB)
      WRITE(*,*) "HD Model: NL =" , NL, "  NB = ", NB
!
!     *** Feld Dimensionierung
      ALLOCATE(MASK_MOU(NL,NB))
      ALLOCATE(FRIV(NL,NB))
      ALLOCATE(INDEXX(NL,NB))
      ALLOCATE(INDEXY(NL,NB))

      ierr = nf90_inq_varid(ncid,'FMOU_HD_TO_NEMO',varid)
      ierr = nf90_get_var(ncid,varid,MASK_MOU)

      ierr = nf90_inq_varid(ncid,'INDEXX',varid)
      ierr = nf90_get_var(ncid,varid,INDEXX)
      ierr = nf90_inq_varid(ncid,'INDEXY',varid)
      ierr = nf90_get_var(ncid,varid,INDEXY)
      ierr = nf90_close(ncid)
      WRITE(*,*) TRIM(DNHDTONEMO)," is closed."
!
      WRITE(*,*) "Anzahl der HD Mouth points", SUM(MASK_MOU, MASK=MASK_MOU.GT.0.5)
!
!!        FLON(JL,JB) = FLOAT(JL)/ 2. - 180.25 
!!        FLAT(JL,JB) = 90.25 - FLOAT(JB) / 2.
!	  
!     ******** READ NEMO grid info 
      ierr = nf90_open(DNNEMO, NF90_NOWRITE, ncid)
      ierr = nf90_inq_dimid(ncid,'x',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
      ierr = nf90_inq_dimid(ncid,'y',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NY)
      ! Testing --> NX = 619, NY=523
      WRITE(*,*) "NEMO:  NX =" , NX, "  NY = ", NY
!
!     *** Feld Dimensionierung
      ALLOCATE(MASK_NEMO(NX,NY))
      ALLOCATE(XNEMO(NX,NY))
      ALLOCATE(YNEMO(NX,NY))
      ALLOCATE(FRIV_NEMO(NX,NY))
!
!     Reading arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'FMOU_HD_ON_NEMO',varid)
      ierr = nf90_get_var(ncid,varid,MASK_NEMO)
      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_get_var(ncid,varid,XNEMO)
      ierr = nf90_inq_varid(ncid,'lat',varid)
      ierr = nf90_get_var(ncid,varid,YNEMO)
      ierr = nf90_close(ncid)
      WRITE(*,*) "Anzahl der NEMO Mouth points", SUM(MASK_NEMO, MASK=MASK_NEMO.GT.0.5)
!
!     ******** Open discharge input and NEMO discarge output files
      ierr = nf90_open(DNINP, NF90_NOWRITE, ncid2)
!
!     ******** Read HD model river discharge
      ierr = nf90_inq_varid(ncid2,'friv',varid)
      ierr = nf90_get_var(ncid2,varid,FRIV)
      ierr = nf90_close(ncid2)
      WHERE (MASK_MOU(:,:).LT.0.5)
        FRIV(:,:) = 0.
      END WHERE
!    
!     ******** Transfer HD model river discharge to NEMO grid
      CALL DIS_TO_NEMO(NL, NB, FRIV, INDEXX, INDEXY, NX, NY, FRIV_NEMO)
      WRITE(*,*) "     Summe HD: ", SUM(FRIV(:,:) * MASK_MOU(:,:))
      WRITE(*,*) "Summe HD-NEMO: ", SUM(FRIV_NEMO(:,:) * MASK_NEMO(:,:))
!
      IF (ILLOG.GT.0) THEN
         WRITE(*,*) 'Log ouput discharge HD:      ', FRIV(ILLOG, IBLOG), ' at ', ILLOG, IBLOG 
         WRITE(*,*) 'Log ouput discharge HD-NEMO: ', FRIV_NEMO(INDEXX(ILLOG, IBLOG),INDEXY(ILLOG, IBLOG)), &
                ' at ', INDEXX(ILLOG, IBLOG),INDEXY(ILLOG, IBLOG)
      ENDIF
!
!     ******** Write river discharge on NEMO grid
!
      DNOUT="dis_nemo.nc"
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'x', NX, dimids(1))
      ierr = nf90_def_dim(ncid, 'y', NY, dimids(2))

      ierr = nf90_def_var(ncid, 'FRIV_NEMO', NF90_float, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 219)
      ierr = nf90_put_att(ncid,varid,'units','[m3 s-1]')
      ierr = nf90_put_att(ncid,varid,'long_name','river discharge')
      ierr = nf90_put_att(ncid,varid,'standard_name','Water_volume_transport_into_NEMO_ocean_from_rivers')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, FRIV_NEMO)

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
      WRITE(*,*) "River-Discharge into NEMO ocean was written: ", TRIM(DNOUT)

      DEALLOCATE(MASK_MOU)
      DEALLOCATE(FRIV)
      DEALLOCATE(INDEXX)
      DEALLOCATE(INDEXY)
      DEALLOCATE(MASK_NEMO)
      DEALLOCATE(FRIV_NEMO)
      DEALLOCATE(XNEMO)
      DEALLOCATE(YNEMO)
!
!     ******** Programmende
  999 END
!
!     *********************************************************************
      SUBROUTINE DIS_TO_NEMO(NL, NB, FRIV, INDEXX, INDEXY, NX, NY, FRIV_NEMO)

!     *********************************************************************
!
!     Transfer HD discharge to NEMO grid

      INTEGER, INTENT(in) :: NL, NB, NX, NY
      REAL, DIMENSION(NL,NB), INTENT(in) :: FRIV
      INTEGER, DIMENSION(NL,NB), INTENT(in) :: INDEXX
      INTEGER, DIMENSION(NL,NB), INTENT(in) :: INDEXY
      REAL, DIMENSION(NX,NY), INTENT(out) :: FRIV_NEMO

      INTEGER :: JL, JB
!
      FRIV_NEMO(:,:) = 0.
      DO JB=1, NB
      DO JL=1, NL 
      IF (INDEXX(JL,JB).GT.0.5) THEN
        FRIV_NEMO(INDEXX(JL,JB), INDEXY(JL,JB)) = FRIV_NEMO(INDEXX(JL,JB), INDEXY(JL,JB)) + FRIV(JL,JB)
      ENDIF
      ENDDO	  
      ENDDO	  
!
      END SUBROUTINE
!	  
