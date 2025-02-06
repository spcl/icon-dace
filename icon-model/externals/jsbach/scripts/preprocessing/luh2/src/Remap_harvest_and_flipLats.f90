!> Program to remap harvested biomass data
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>

! Program taken from the script collection of SW/TR (rewrite_lu_harvest.f90) to remap
! harvested biomass data which was previously aggregated from LUH2_v1 data from 0.25deg to ires
!
!-- JN 01.07.16:
!* changed the expected input grid
!* new files are in kgC not in MgC
!* added reading of float variables
!-- TR 25.08.16:
!* processing of years with 3 digits
!-- JN 06.12.16:
!* also read work path from timeinfo.nml
!
PROGRAM remap_harvest_and_flip
   IMPLICIT NONE
   INCLUDE 'netcdf.inc'
   ! Number model from which the SELECTED_*_KIND are requested:
   !
   !                   4 byte REAL      8 byte REAL
   !          CRAY:        -            precision =   13
   !                                    exponent  = 2465
   !          IEEE:    precision =  6   precision =   15
   !                   exponent  = 37   exponent  =  307
   INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
   INTEGER, PARAMETER :: wp = dp   ! working precision
   INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
   INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
! error status return
   INTEGER  ::  IRET
! netCDF id
   INTEGER  ::  NCID,NCIDOUT,NVARID
! dimensions
   INTEGER  ::  NDIMLON,NDIMLAT,NDIMID
   INTEGER  ::  NDIM,NVAR,NATTR,NUNLIMDIM
   INTEGER  ::  NDIMINID(3), NDIMOUTID(3), NVAROUTID(4)
   INTEGER, ALLOCATABLE, DIMENSION(:)  ::  NDIMLEN
   CHARACTER(LEN=60), ALLOCATABLE, DIMENSION(:)  ::  FDIMNAM
! variables
   INTEGER  ::  NIN(100)
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARTYP
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARDIM
   INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: NVARDIMID
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARATT
   CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:)  ::  FVARNAM
   REAL(dp), ALLOCATABLE, DIMENSION(:)  :: RVARTIME
! attributes
   INTEGER  :: NVARATTMAX
   CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:,:)  ::  FATTNAM
   CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:,:)  ::  FATT
   INTEGER, ALLOCATABLE, DIMENSION(:,:)  ::  NATTTYP
   INTEGER, ALLOCATABLE, DIMENSION(:,:)  ::  NATTDIM
   INTEGER, ALLOCATABLE, DIMENSION(:,:)  ::  NATT
   REAL(dp), ALLOCATABLE, DIMENSION(:,:)  ::  RATT
! data
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RDATA
! gauss output
   INTEGER  ::  NGLON, NGLAT, ISTATES_OUT
   REAL(dp)  ::  UMFANG
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  ROUTPUT_GAUSS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RTRANS_OUT_GAUSS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:)  ::  RTRANS_FLIP
   REAL(dp),ALLOCATABLE, DIMENSION(:,:)  ::  RSTATES_FLIP
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XILON
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XILAT
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLON
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLAT
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLAT_SWITCH
! indices
   INTEGER  ::  I,II,J,ITIME
   INTEGER  :: year_start, year_end, nyear, ires
! work path
   CHARACTER(LEN=128) :: workingPath
! file names
   CHARACTER(LEN=128), parameter  ::  INPUT_FILES= "/lu_out_harvest/" ! directory of input data files for land use change
   CHARACTER(LEN=128), parameter  ::  OUTPUT_FILES= "/lu_out_harvest/" ! directory of output data files
!!$   CHARACTER(LEN=128), parameter  ::  INPUT_FILES= "/future_image/lu_out_harvest/" ! directory of input data files
!!$   CHARACTER(LEN=128), parameter  ::  OUTPUT_FILES= "/future_image/lu_out_harvest/" ! directory of output data files
   CHARACTER(LEN=128)  ::  FNAM
   CHARACTER(LEN=4)  ::  INPUT_YEAR
   CHARACTER(LEN=3)  ::  FRES
! parameters
   INTEGER, PARAMETER  ::  NLON = 1440            ! number of longitudes
   INTEGER, PARAMETER  ::  NLAT = 720            ! number of latitudes
! logicals
   LOGICAL  ::  LLON, LLAT
   LOGICAL  ::  LREAD, LWRITE
   LOGICAL  ::  INFO = .FALSE.

   NAMELIST /timeinfo/ nyear, year_start, year_end, ires, workingPath

   ! Read time info from namelist
   ires       = 63
   nyear      = -1
   year_start = -1
   year_end   = -2
   workingPath = "."
   open(11,file='timeinfo.nml',status='old',position='rewind')
   read(11,timeinfo,iostat=i)
   if (i /= 0) then
      write(*,*) 'Namelist not found or corrupted'
      stop
   endif
   close(11)
   if ((year_end < year_start) .and. (nyear < 0)) then
     write(*,*) 'Bad combination of namelist parameters'
     stop
   endif
   if (year_end < year_start) year_end = year_start + nyear - 1

   year_start = year_start - 1            ! To be adjusted in the time loop
   nyear      = year_end - year_start

   write(*,*) 'Work path: ',trim(workingPath)
   write(*,*) 'Processing years: ',year_start+1,year_end

   DO ITIME = 1,NYEAR  ! open time loop

      IF (year_start + ITIME .LE. 999) THEN
         WRITE(INPUT_YEAR,'(A,i3)') '0' , year_start + ITIME
      ELSE
         WRITE(INPUT_YEAR,'(i4)') year_start + ITIME
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! INPUT-FILE: Harvest !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      LREAD=.FALSE.
      FNAM = trim(workingPath) // trim(INPUT_FILES) // 'LUH_harvest_' // TRIM(INPUT_YEAR) // '.nc'

! open the input-file
      WRITE (*,*) 'FNAM: ',TRIM(FNAM)
      IRET = nf_open(FNAM,nf_nowrite,NCID)
      CALL check_err(IRET)
! check what is in the input-file
      IRET = nf_inq(NCID,NDIM,NVAR,NATTR,NUNLIMDIM)
      CALL check_err(IRET)
      IF (INFO) WRITE (*,*) 'NDIM,NVAR,NATTR,NUNLIMDIM',NDIM,NVAR,NATTR,NUNLIMDIM
! get the dimension name and length
      ALLOCATE (NDIMLEN(NDIM))
      ALLOCATE (FDIMNAM(NDIM))
      DO I=1,NDIM
        IRET = nf_inq_dim(NCID,I,FDIMNAM(I),NDIMLEN(I))
        CALL check_err(IRET)
        IF (INFO) WRITE (*,*) 'I,FDIMNAM,NDIMLEN',I,TRIM(FDIMNAM(I)),NDIMLEN(I)
      END DO
! get variable names, types and shapes
      ALLOCATE (FVARNAM(NVAR))
      ALLOCATE (NVARTYP(NVAR))
      ALLOCATE (NVARDIM(NVAR))
      ALLOCATE (NVARDIMID(NVAR,100))
      ALLOCATE (NVARATT(NVAR))
      NIN(:) = 0
      DO I=1,NVAR
       IRET = nf_inq_var(NCID,I,FVARNAM(I),NVARTYP(I),NVARDIM(I),NIN,NVARATT(I))
        CALL check_err(IRET)
        DO II=1,100
          NVARDIMID(I,II)=NIN(II)
        END DO
        IF (INFO) WRITE (*,*) 'I,FVARNAM,NVARTYP,NVARDIM,NVARDIMID,NVARATT', &
           I,TRIM(FVARNAM(I)),NVARTYP(I),NVARDIM(I),NVARDIMID(I,1),NVARATT(I)
      END DO
! get attribute names and then attribute types and lengths
      NVARATTMAX=0
      DO I=1,NVAR
        IF (NVARATT(I).GT.NVARATTMAX) NVARATTMAX=NVARATT(I)
      END DO
      IF (INFO) WRITE (*,*) 'NVARATTMAX', NVARATTMAX
      ALLOCATE (FATTNAM(NVAR,NVARATTMAX))
      ALLOCATE (NATTTYP(NVAR,NVARATTMAX))
      ALLOCATE (NATTDIM(NVAR,NVARATTMAX))
      DO I=1,NVAR
        DO II=1,NVARATT(I)
          IRET = nf_inq_attname(NCID,I,II,FATTNAM(I,II))
          CALL check_err(IRET)
          IRET = nf_inq_att(NCID,I,FATTNAM(I,II),NATTTYP(I,II),NATTDIM(I,II))
          CALL check_err(IRET)
          IF (INFO) WRITE (*,*) 'I,II,FATTNAME,ATTTYP,NATTDIM',I,II,TRIM(FATTNAM(I,II)),NATTTYP(I,II),NATTDIM(I,II)
        END DO
      END DO
! get attributes
      ALLOCATE (FATT(NVAR,NVARATTMAX))
      ALLOCATE (NATT(NVAR,NVARATTMAX))
      ALLOCATE (RATT(NVAR,NVARATTMAX))
      FATT(:,:) = ''
      DO I=1,NVAR
        DO II=1,NVARATT(I)
          IF (NATTTYP(I,II).EQ.2) THEN
            IRET = nf_get_att_text(NCID,I,FATTNAM(I,II),FATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,FATT',I,TRIM(FATTNAM(I,II)),TRIM(FATT(I,II))
            IF (TRIM(FVARNAM(I))=='harvest' .AND. TRIM(FATTNAM(I,II))=='units') THEN
               NATTDIM(I,II)=12
               FATT(I,II)='mol(C)m-2s-1'
            END IF
          ELSE IF (NATTTYP(I,II).EQ.4) THEN
            IRET = nf_get_att_int(NCID,I,FATTNAM(I,II),NATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,NATT',I,TRIM(FATTNAM(I,II)),NATT(I,II)
          ELSE IF ((NATTTYP(I,II).EQ.5) .OR. (NATTTYP(I,II).EQ.6)) THEN
            IRET = nf_get_att_double(NCID,I,FATTNAM(I,II),RATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,RATT',I,TRIM(FATTNAM(I,II)),RATT(I,II)
          ELSE
            WRITE (*,*) 'Type of attribute not prepared in the code, failed to read grid',I,II,NATTTYP(I,II)
            STOP
          END IF
        END DO
      END DO
! get data
      DO I = 1,NVAR
         IF (.NOT. FVARNAM(I) == 'lon' .AND. .NOT. FVARNAM(I) == 'lat' .AND. .NOT. FVARNAM(I) == 'time') THEN
            IF (LREAD) STOP 'More data fields in the file than expected!'
            LREAD=.TRUE.
            ALLOCATE(RDATA(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),1))
            IRET = nf_get_var_double(NCID,I,RDATA)
            CALL check_err(IRET)
         ELSE IF (FVARNAM(I) == 'time') THEN
            ALLOCATE(RVARTIME(NDIMLEN(NVARDIMID(I,1))))
            IRET = nf_get_var_double(NCID,I,RVARTIME)
            CALL check_err(IRET)
         END IF
      END DO

      IF (.NOT. ALLOCATED(RDATA)) STOP 'No data field found!'
      If (.NOT. ALLOCATED(RVARTIME)) STOP 'No time step found!'

      WHERE (RDATA(:,:,1) > 1.e+12)
         RDATA(:,:,1) = -9999.
      END WHERE

! close the input-file states
      IRET = nf_close(NCID)
      CALL check_err(IRET)

! check number of longitudes and latitudes
      LLON = .false.
      LLAT = .false.
      DO I = 1,NDIM
         IF (FDIMNAM(I) == 'lon') THEN
            NDIMINID(1) = I
            IF (NDIMLEN(I) == NLON) LLON = .true.
         ELSEIF (FDIMNAM(I) == 'lat') THEN
            NDIMINID(2) = I
            IF (NDIMLEN(I) == NLAT) LLAT = .true.
         ELSE IF (FDIMNAM(I) == 'time') THEN
            NDIMINID(3) = I
         END IF
      END DO

      IF (.NOT. LLON) STOP 'Dimension lon not found or number of longitudes not correct'
      IF (.NOT. LLAT) STOP 'Dimension lat not found or number of latitudes not correct'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! OUTPUT TRANSITIONS GAUSS GRID !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IRES <= 99) THEN
         WRITE(FRES,'(i2)')  IRES
      ELSE IF (IRES >=100 .AND. IRES <= 999) THEN
         WRITE(FRES,'(i3)')  IRES
      ELSE
         STOP 'Resolution not prepared!'
      END IF

      IF (INFO) WRITE (*,*) 'FRES,IRES',FRES,IRES
      IF (IRES == 21) NGLAT = 32
      IF (IRES == 31) NGLAT = 48
      IF (IRES == 42) NGLAT = 64
      IF (IRES == 63) NGLAT = 96
      IF (IRES == 85) NGLAT = 128
      IF (IRES == 106) NGLAT = 160
      IF (IRES == 127) NGLAT = 192
      IF (IRES == 159) NGLAT = 240
      IF (IRES == 255) NGLAT = 384
      IF (IRES == 319) NGLAT = 480
      NGLON = NGLAT * 2
      UMFANG = 360.

!  calculate incoming latitude
      IF (INFO) WRITE(*,*) 'INCOMING GRID:'
      ALLOCATE(XILAT(NLAT+2))
      ALLOCATE(XILON(NLON+2))
      DO J = 1,NLAT+2
         XILAT(J) = - UMFANG / (4. * FLOAT(NLAT)) - UMFANG / 4. + UMFANG / 2. * FLOAT(J-1) / FLOAT(NLAT)
      END DO
      IF (INFO) WRITE(*,'(A,3F8.3,A,2F8.3)') ' XILAT : ',XILAT(1),XILAT(2),XILAT(3),' ... ',XILAT(NLAT+1),XILAT(NLAT+2)
!  calculate incoming longitude
      DO I = 1,NLON+2
         XILON(I) = (FLOAT(I-1) - 0.5) * UMFANG / FLOAT(NLON)
      END DO
      IF (INFO) WRITE(*,'(A,3F8.3,A,2F8.3)') ' XILON : ',XILON(1),XILON(2),XILON(3),' ... ',XILON(NLON+1),XILON(NLON+2)
!  calculate outgoing latitude
      IF (INFO) WRITE(*,*) 'OUTGOING GRID:'
      ALLOCATE (XGLAT(NGLAT))
      ALLOCATE (XGLAT_SWITCH(NGLAT))
      ALLOCATE (XGLON(NGLON))
      CALL GAUAW (XGLAT,NGLAT)
      DO J = 1,NGLAT
         XGLAT(J) = -ASIN(XGLAT(J)) / (8. * ATAN(1.)) * UMFANG
      END DO
      IF (INFO) WRITE(*,'(A,3F8.3,A,2F8.3)') ' XGLAT  : ',XGLAT(1),XGLAT(2),XGLAT(3),' ... ',XGLAT(NGLAT-1),XGLAT(NGLAT)
!  calculate outgoing longitude
      DO I = 1,NGLON
         XGLON(I) = (I - 1.) / FLOAT(NGLON) * UMFANG
      END DO
      IF (INFO) WRITE(*,'(A,3F8.3,A,2F8.3)') ' XGLON  : ',XGLON(1),XGLON(2),XGLON(3),' ... ',XGLON(NGLON-1),XGLON(NGLON)

      ALLOCATE(RTRANS_FLIP(NLON+2,NLAT+2))
      ALLOCATE(RTRANS_OUT_GAUSS(NGLON,NGLAT,1))
      ALLOCATE(ROUTPUT_GAUSS(NGLON,NGLAT,1))

      RTRANS_OUT_GAUSS(:,:,:) = 0.
      RTRANS_FLIP(:,:) = 0.
      DO J = 1,NLAT
         DO I = 1,NLON/2
             RTRANS_FLIP(I + 1,J + 1) = RDATA(I + NLON/2,J,1)
         END DO
         DO I = NLON/2 + 1,NLON
            RTRANS_FLIP(I + 1,J + 1) = RDATA(I - NLON/2,J,1)
         END DO
         RTRANS_FLIP(1,J + 1) = RTRANS_FLIP(NLON + 1,J + 1)
         RTRANS_FLIP(NLON + 2, J + 1) = RTRANS_FLIP(2,J + 1)
      END DO
      DO I = 1,NLON + 2
         RTRANS_FLIP(I,1) = RTRANS_FLIP(I,2)
         RTRANS_FLIP(I,NLAT + 2) = RTRANS_FLIP(I,NLAT + 1)
      END DO
      IF (INFO) WRITE (*,*) 'call intpol2'
      CALL INTPOL2(INFO,NLON+2,NLAT+2,RTRANS_FLIP,XILON,XILAT,NGLON,NGLAT,XGLON,XGLAT,RTRANS_OUT_GAUSS(:,:,1))
!!$   FNAM = trim(workingPath) // trim(OUTPUT_FILES) // 'LUH_harvest_T' // TRIM(FRES) // '_rcp26_' // TRIM(INPUT_YEAR) // '.nc'
      FNAM = trim(workingPath) // trim(OUTPUT_FILES) // 'LUH_harvest_T' // TRIM(FRES) // '_' // TRIM(INPUT_YEAR) // '.nc'
! switch latitudes
      DO J = 1,NGLAT
         XGLAT_SWITCH(J) = XGLAT(NGLAT + 1 -J)
         ROUTPUT_GAUSS(:,J,1) = RTRANS_OUT_GAUSS(:,NGLAT+1-J,1)
      END DO
! open the output-file
      IF (INFO) WRITE (*,*) 'Output to: ',TRIM(FNAM)
      IRET = nf_create(FNAM,nf_clobber,NCIDOUT)
      CALL check_err(IRET)
! define dimensions (same as the dimensions in the transitions input-file)
      IF (INFO) WRITE (*,*) 'Define dimensions: ',NGLON,NGLAT,NDIM
      NDIMLEN(:) = (/NGLON,NGLAT,NF_UNLIMITED/) ! lon, lat, time
      DO I = 1,NDIM
         CALL check_err(nf_def_dim(NCIDOUT,FDIMNAM(NDIMINID(I)),NDIMLEN(I),NDIMOUTID(I)))
         CALL check_err(nf_def_var(NCIDOUT,FVARNAM(NDIMINID(I)),6,1,NDIMOUTID(I),NVAROUTID(I)))
         CALL PUTATT(NVAR,NDIMINID(I),NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(I),NATT,RATT,INFO)
      END DO
      CALL check_err(nf_def_var(NCIDOUT,'harvest',6,3,NDIMOUTID,NVAROUTID(4)))
      CALL PUTATT(NVAR,4,NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(4),NATT,RATT,INFO)
      CALL check_err(nf_enddef(NCIDOUT))
      CALL check_err(IRET)
      IF (INFO) WRITE (*,*) 'Dimensions and variables defined.'
! Write data
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(1),XGLON))
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(2),XGLAT_SWITCH))
      CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(3),(/1/),(/1/),RVARTIME))
      CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(4),(/1,1,1/),(/NGLON,NGLAT,1/),ROUTPUT_GAUSS(:,:,1)))
      IF (INFO) WRITE (*,*) 'Data wrote.'

! close the output-file
      CALL check_err(nf_close(NCIDOUT))

! deallocate arrays
      DEALLOCATE(NDIMLEN,FDIMNAM,FVARNAM,NVARTYP,NVARDIM,NVARDIMID,NVARATT)
      DEALLOCATE(FATTNAM,NATTTYP,NATTDIM,FATT,NATT,RATT)
      DEALLOCATE(RDATA,RVARTIME)
      DEALLOCATE(XILAT,XILON,XGLAT,XGLON,XGLAT_SWITCH)
      DEALLOCATE(RTRANS_FLIP,RTRANS_OUT_GAUSS,ROUTPUT_GAUSS)

   END DO  ! close time loop

   WRITE (*,*) "Program sucessfully finished."

CONTAINS

  ! --- check_err --------------------------------------------------------------------------------------------

   SUBROUTINE check_err(iret,text)
   IMPLICIT NONE
   INCLUDE 'netcdf.inc'
   INTEGER,INTENT(IN)          :: iret
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: text

      IF (iret /= NF_NOERR) THEN
         IF (present(text)) WRITE(*,*) 'check_err(): '//TRIM(text)
         WRITE(*,*) 'check_err(): netcdf says: '//TRIM(nf_strerror(iret))
         STOP
      END IF

   END SUBROUTINE check_err

   SUBROUTINE PUTATT(NVAR,I,NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVARID,NATT,RATT,INFO)
   !  This routine adds the attributes to the output-file
      IMPLICIT NONE
      INCLUDE 'netcdf.inc'
      INTEGER, INTENT(IN)  ::  NVAR
      INTEGER, INTENT(IN)  ::  I
      INTEGER, INTENT(IN), DIMENSION(NVAR)  ::  NVARATT
      INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I))  :: NATTTYP
      CHARACTER(LEN=*), INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  FATTNAM
      INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I)) ::  NATTDIM
      CHARACTER(LEN=*), INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  FATT
      INTEGER, INTENT(IN)  ::  NCIDOUT, NVARID
      INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  NATT
      REAL(dp), INTENT(IN), DIMENSION(NVAR,NVARATT(I)) ::  RATT
      LOGICAL, INTENT(IN)  ::  INFO
      INTEGER  ::  III,IRET

      DO III=1,NVARATT(I)
        IF (NATTTYP(I,III).EQ.2) THEN
          IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTDIM(I,III),TRIM(FATT(I,III))
          IRET = nf_put_att_text(NCIDOUT,NVARID,FATTNAM(I,III),NATTDIM(I,III),FATT(I,III))
          CALL check_err(IRET)
        ELSE IF (NATTTYP(I,III).EQ.4) THEN
           IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTTYP(I,III),NATTDIM(I,III),NATT(I,III)
           IRET = nf_put_att_int(NCIDOUT,NVARID,FATTNAM(I,III),NATTTYP(I,III),NATTDIM(I,III),NATT(I,III))
           CALL check_err(IRET)
        ELSE IF ((NATTTYP(I,III).EQ.5) .OR. (NATTTYP(I,III).EQ.6)) THEN
           IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTTYP(I,III),NATTDIM(I,III),RATT(I,III)
           IRET = nf_put_att_double(NCIDOUT,NVARID,FATTNAM(I,III),NATTTYP(I,III),NATTDIM(I,III),RATT(I,III))
           CALL check_err(IRET)
         ELSE
           WRITE (*,*) 'Type of attribute not prepared in the code, failed to write',I,III,NATTTYP(I,III)
         END IF
      END DO

   END SUBROUTINE PUTATT


SUBROUTINE INTPOL2(INFO,NXM,NYM,FIELDM,XM,YM,NX,NY,XX,YY,FIELD)
! Flux conserving transformation of a global field with a high resolution to a lower resolution grid.
! The incoming field is assumed to range from the South Pole to the North Pole.
! This routine is designed for pre- and post-processing purposes It is not optimized in terms of computational efficiency.
! ATTENTION: If the outgoing grid is finer than the incoming grid the result will be erroneous without warning!
! author: Thomas Raddatz 2003 (at MPI-MET Hamburg)
   IMPLICIT NONE
   LOGICAL,INTENT(IN) :: INFO
   INTEGER,INTENT(IN) :: NXM,NYM
   INTEGER,INTENT(IN) :: NX,NY
   REAL(dp),INTENT(IN)    :: FIELDM(NXM,NYM)
   REAL(dp),INTENT(IN)    :: XM(NXM),YM(NYM)
   REAL(dp),INTENT(IN)    :: XX(NX),YY(NY)
   REAL(dp),INTENT(OUT)   :: FIELD(NX,NY)
! local variables
   REAL(dp),PARAMETER :: deg_to_rad=3.1415926536/180.
   REAL(dp),PARAMETER :: EARTH_RADIUS=6371001.
   REAL(dp),PARAMETER :: SECONDS_IN_YEAR=31557600. ! this is 365.25 days (including leap years)
   REAL(dp),PARAMETER :: kgC_in_molC=1.e+03/12.011
   REAL(dp) :: FIELD_FLUX
   REAL(dp) :: AREA(NX,NY),XDEL,AIN,AIN1,AIN2,AIN3,AIN4,PI
   REAL(dp) :: X(NX+3),Y(NY+2)
   INTEGER :: IFLAG,IX,IY,ICASE
! indices
   INTEGER :: I,J,II,JJ
   PI = 3.14159265
! create extended arrays x,y
   DO I=2,NX+1
     X(I)=XX(I-1)
   END DO
   DO J=2,NY+1
     Y(J) = YY(J-1)
   END DO
   X(1) = (2. * X(2)) - X(3)
   X(NX+2) = (2. * X(NX+1)) - X(NX)
   X(NX+3) = (2. * X(NX+2)) - X(NX+1)
   Y(1) = -180. - Y(2) - 1.E-04
   Y(NY+2) = 180. - Y(NY+1) + 1.E-04
   IF (INFO) THEN
     WRITE(*,*) 'INCOMING GRID:'
     WRITE(*,'(A,3F8.3,A,2F8.3)') ' YM : ',YM(1),YM(2),YM(3),' ... ',YM(NYM-1),YM(NYM)
     WRITE(*,'(A,3F8.3,A,2F8.3)') ' XM : ',XM(1),XM(2),XM(3),' ... ',XM(NXM-1),XM(NXM)
     WRITE(*,*) 'OUTGOING GRID:'
     WRITE(*,'(A,3F8.3,A,2F8.3)') ' Y  : ',Y(1),Y(2),Y(3),' ... ',Y(NY+1),Y(NY+2)
     WRITE(*,'(A,3F8.3,A,2F8.3)') ' X  : ',X(1),X(2),X(3),' ... ',X(NX+2),X(NX+3)
   END IF
! initialise FIELD and AREA
   FIELD(:,:) = 0.
   AREA(:,:) = 0.
! loop: All IN-GRID POINTS
   XDEL = ABS(XM(1) - XM(2))
   DO 99 II=2,NXM-1
     DO 98 JJ=2,NYM-1
       IF (FIELDM(II,JJ) > -1.) THEN   ! exclude ocean grid points
       IFLAG=0
       AIN = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.)) * XDEL
       FIELD_FLUX = FIELDM(II,JJ) * kgC_in_molC / & ! translate absolute value [MgC] in flux [mol(C)m-2s-1]
                    (AIN * EARTH_RADIUS * EARTH_RADIUS * SECONDS_IN_YEAR * PI / 180.)
! assign lower left corner of the IN GRID BOX to one OUT GRID BOX
       DO I=2,NX+2
         DO J=2,NY+1
           IF (XM(II) + XM(II-1) >= X(I) + X(I-1) .AND. XM(II) + XM(II-1) < X(I) + X(I+1) .AND. &
               YM(JJ) + YM(JJ-1) >= Y(J) + Y(J-1) .AND. YM(JJ) + YM(JJ-1) < Y(J) + Y(J+1)) THEN
             IFLAG = IFLAG + 1
             IX = I
             IY = J
	   ENDIF
         ENDDO
       ENDDO
       IF (IFLAG == 0) STOP 'FAILED TO FIND IX,IY'
       IF (IFLAG >= 2) STOP 'FOUND MULTIPLE IX,IY'
! 4 cases:
       IFLAG = 0
       IF (XM(II) + XM(II+1) < X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) < Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX is completely in ONE OUT-GRID BOX
         IFLAG = IFLAG + 1
         ICASE = 1
       ELSEIF (XM(II) + XM(II+1) < X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) >= Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX stretches over 2 OUT-GRID BOXES in meridional direction
         IFLAG = IFLAG + 1
         ICASE = 2
       ELSEIF (XM(II) + XM(II+1) >= X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) < Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX stretches over 2 OUT-GRID BOXES in longitudinal direction
         IFLAG = IFLAG+1
         ICASE = 3
       ELSEIF (XM(II) + XM(II+1) >= X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) >= Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX partially covers 4 OUT-GRID BOXES
         IFLAG = IFLAG + 1
         ICASE = 4
       ENDIF
       IF (IX == NX+2) IX = 2
       IF (IFLAG == 0) STOP 'FAILED TO FIND ICASE'
       IF (IFLAG.GE.2) STOP 'FOUND MULTIPLE ICASE'
! add flux multiplied by the associated area to OUT-GRID
       IF (ICASE == 1) THEN
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + AIN
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELD_FLUX * AIN)
       ELSEIF (ICASE == 2) THEN
         IF (IY > NY) WRITE(*,*) 'IY,NY: ATTENTION (ICASE=2)!',IY,NY
         AIN1 = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * XDEL
         AIN2 = ABS(SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * XDEL
         IF (ABS(AIN1 + AIN2 - AIN) / AIN > 1.E-04) STOP 'AREAS 2 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + AIN1
         AREA(IX-1,IY) = AREA(IX-1,IY) + AIN2
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELD_FLUX * AIN1)
         FIELD(IX-1,IY) = FIELD(IX-1,IY) + (FIELD_FLUX * AIN2)
       ELSEIF (ICASE == 3) THEN
         AIN1 = AIN * 0.5 * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / XDEL
         AIN2 = AIN * 0.5 * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / XDEL
         IF (ABS(AIN1 + AIN2 - AIN) / AIN > 1.E-04) STOP 'AREAS 3 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + AIN1
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELD_FLUX * AIN1)
         IF (IX == NX+1) THEN
           AREA(1,IY-1) = AREA(1,IY-1) + AIN2
           FIELD(1,IY-1) = FIELD(1,IY-1) + (FIELD_FLUX * AIN2)
         ELSE
           AREA(IX,IY-1) = AREA(IX,IY-1) + AIN2
           FIELD(IX,IY-1) = FIELD(IX,IY-1) + (FIELD_FLUX * AIN2)
         ENDIF
       ELSEIF (ICASE.EQ.4) THEN
         IF (IY > NY) WRITE (*,*) 'IY,NY: ACHTUNG (ICASE=4)!',IY,NY
!!$         AIN1 = ABS(SIND((YM(JJ-1) + YM(JJ)) / 2.) - SIND((Y(IY) + Y(IY+1)) / 2.)) * &
!!$                XDEL * 90. * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * PI)
!!$         AIN2 = ABS(SIND((YM(JJ-1) + YM(JJ)) / 2.) - SIND((Y(IY) + Y(IY+1)) / 2.)) * &
!!$                XDEL * 90. * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * PI)
!!$         AIN3 = ABS(SIND((YM(JJ) + YM(JJ+1)) / 2.) - SIND((Y(IY) + Y(IY+1)) / 2.)) * &
!!$                XDEL * 90. * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * PI)
!!$         AIN4 = ABS(SIND((YM(JJ) + YM(JJ+1)) / 2.) - SIND((Y(IY) + Y(IY+1)) / 2.)) * &
!!$                XDEL * 90. * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * PI)

         AIN1 = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * 2.)
         AIN2 = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * 2.)
         AIN3 = ABS(SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * 2.)
         AIN4 = ABS(SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * 2.)

         IF (ABS(AIN1 + AIN2 + AIN3 + AIN4 - AIN) / AIN > 1.E-04) STOP 'AREAS 4 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + AIN1
         AREA(IX-1,IY) = AREA(IX-1,IY) + AIN3
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELD_FLUX * AIN1)
         FIELD(IX-1,IY) = FIELD(IX-1,IY) + (FIELD_FLUX * AIN3)
         IF (IX == NX+1) THEN
           AREA(1,IY-1) = AREA(1,IY-1) + AIN2
           AREA(1,IY) = AREA(1,IY) + AIN4
           FIELD(1,IY-1) = FIELD(1,IY-1) + (FIELD_FLUX * AIN2)
           FIELD(1,IY) = FIELD(1,IY) + (FIELD_FLUX * AIN4)
         ELSE
           AREA(IX,IY-1) = AREA(IX,IY-1) + AIN2
           AREA(IX,IY) = AREA(IX,IY) + AIN4
           FIELD(IX,IY-1) = FIELD(IX,IY-1) + (FIELD_FLUX * AIN2)
           FIELD(IX,IY) = FIELD(IX,IY) + (FIELD_FLUX * AIN4)
         ENDIF
       ENDIF
     END IF
98   CONTINUE
99 CONTINUE
! divide by area
   WHERE (AREA(:,:) > 1.e-12)
     FIELD(:,:) = FIELD(:,:) / AREA(:,:)
   ELSEWHERE
     FIELD(:,:) = 0. ! FILL_VALUE
   END WHERE
END SUBROUTINE INTPOL2

SUBROUTINE GAUAW (pa,nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    IMPLICIT NONE

    !  Scalar arguments
    INTEGER :: nlat

    !  Array arguments
    REAL(dp)    :: pa(nlat)
    ! *pa*  - array, length at least *k,* to receive abscis abscissas.

    !  Local scalars:
    REAL, PARAMETER :: eps = EPSILON(0.0)
    INTEGER, PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL    :: api
    REAL    :: za, zw, z, zan
    REAL    :: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions
    INTRINSIC ABS, COS, MOD, TAN

    !  Executable statements

    api    = 2.*ASIN(1.)

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat

    DO jgl = 1, ins2
       z = REAL(4*jgl-1)*api/REAL(4*nlat+2)
       pa(jgl) = COS(z+1./(TAN(z)*REAL(8*nlat**2)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)

       DO iter = 1, itemax+1
          zk = 0.0

          ! Newton iteration step

          zkm2 = 1.0
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1)*zx*zkm1-REAL(jn-1)*zkm2)/REAL(jn)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat)*(zkm1-zx*zk))/(1.-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn

          ! computes weight

          zkm2 = 1.0
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1)*zx*zkm1-REAL(jn-1)*zkm2)/REAL(jn)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0-zx*zx)/(REAL(nlat*nlat)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= eps) EXIT
       END DO
       pa(jgl) = zan
    ENDDO
    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
    ENDDO


END SUBROUTINE GAUAW

END PROGRAM remap_harvest_and_flip
