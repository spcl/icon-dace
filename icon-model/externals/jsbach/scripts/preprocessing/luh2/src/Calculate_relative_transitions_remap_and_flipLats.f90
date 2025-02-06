!> Program to calculate relative from absolute transition data
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

! Program taken from the script collection of SW/TR (rewrite_transfer_matrix.f90) to
! calculate relative from absolute transition data which was previously aggregated from LUH2_v1 data
! and then remap it from 0.25 degC to ires
!
!-- JN 07.07.16:
!* changed the expected input grid
!* removed non existing transitions
!* added reading of floats (write as double)
!-- TR 25.08.16:
!* processing of years with 3 digits
!-- JN 06.12.16:
!* also read work path from timeinfo.nml
!-- JN 22-02-17
!* started to speed up the program
PROGRAM Calculate_relative_transitions_remap_and_flipLats
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
! parameters
   INTEGER, PARAMETER  ::  NTRANS = 9           ! number of transitions
   INTEGER, PARAMETER  ::  NTRANS_OUT = 6        ! number of transitions of the new transition matrix
   INTEGER, PARAMETER  ::  NSTATES = 4           ! number of states
   INTEGER, PARAMETER  ::  NLON = 1440           ! number of longitudes
   INTEGER, PARAMETER  ::  NLAT = 720            ! number of latitudes
   REAL(dp), PARAMETER  ::  MISS_VALUE = NF_FILL_DOUBLE  ! missing value
! error status return
   INTEGER  ::  IRET
! netCDF id
   INTEGER  ::  NCID,NCIDOUT,NVARID
! dimensions
   INTEGER  ::  NDIMLON,NDIMLAT,NDIMID
   INTEGER  ::  NDIM,NVAR,NATTR,NUNLIMDIM
   INTEGER  ::  NDIMINID(3), NDIMOUTID(3), NVAROUTID(3+NTRANS_OUT)
   INTEGER, ALLOCATABLE, DIMENSION(:)  ::  NDIMLEN
   CHARACTER(LEN=60), ALLOCATABLE, DIMENSION(:)  ::  FDIMNAM
! variables
   INTEGER  ::  NIN(100)
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARTYP
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARDIM
   INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: NVARDIMID
   INTEGER, ALLOCATABLE, DIMENSION(:)  :: NVARATT
   CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:)  ::  FVARNAM
   REAL(dp), ALLOCATABLE, DIMENSION(:)  :: RVARLON
   REAL(dp), ALLOCATABLE, DIMENSION(:)  :: RVARLAT
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
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RTRANS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RTRANS_REL
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RTRANS_OUT
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  ROUTPUT
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RSTATES
! no pasture
   REAL(dp)  :: RNAT_TO_NAT, RSTATE_NATURAL_NEXT, RSTATE_MANAGED_NOW
   REAL(dp)  :: RTRANS_NATURAL_MANAGED, RTRANS_MANAGED_NATURAL
! gauss output
   INTEGER  ::  NGLON, NGLAT, ISTATES_OUT
   REAL(dp)  ::  UMFANG, EXCESS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  ROUTPUT_GAUSS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:,:)  ::  RTRANS_OUT_GAUSS
   REAL(dp),ALLOCATABLE, DIMENSION(:,:)  ::  RTRANS_FLIP
   REAL(dp),ALLOCATABLE, DIMENSION(:,:)  ::  RSTATES_FLIP
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XILON
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XILAT
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLON
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLAT
   REAL(dp),ALLOCATABLE, DIMENSION(:)  ::  XGLAT_SWITCH
   REAL(dp),POINTER, DIMENSION(:,:) :: MAX_EXCESS
! indices
   INTEGER  ::  I,II,J,K,ITIME,ITRANS,ISTATES
   INTEGER  ::  ICROP2PASTURE, IPASTURE2CROP, IPRIMARY2PASTURE
   INTEGER  ::  IPRIMARY2CROP, ISECONDARY2CROP, ICROP2SECONDARY
   INTEGER  ::  ISECONDARY2PASTURE, IPASTURE2SECONDARY
   INTEGER  ::  IPRIMARY,ISECONDARY,ICROP,IPASTURE
   INTEGER  ::  year_start, year_end, nyear, ires
! work path
   CHARACTER(LEN=128) :: workingPath
! file names
   CHARACTER(LEN=128), parameter  ::  INPUT_FILES1= "/updated_states_out/" ! directory of input data files for updated states
   CHARACTER(LEN=128), parameter  ::  INPUT_FILES2= "/lu_out/" ! directory of input data files for land use change
   CHARACTER(LEN=128), parameter  ::  OUTPUT_FILES= "/lu_out_gauss/" ! directory of output data files
   CHARACTER(LEN=128)  ::  FNAM
   CHARACTER(LEN=4)  ::  INPUT_YEAR
   CHARACTER(LEN=3)  ::  FRES
! other names
   CHARACTER(LEN=128),ALLOCATABLE, DIMENSION(:)  ::  FTRANSNAMES
   CHARACTER(LEN=128),ALLOCATABLE, DIMENSION(:)  ::  FSTATENAMES
! logicals
   LOGICAL  ::  LLON, LLAT
   LOGICAL  ::  LCROP2PASTURE, LPASTURE2CROP, LPRIMARY2PASTURE
   LOGICAL  ::  LPRIMARY2CROP, LSECONDARY2CROP, LCROP2SECONDARY
   LOGICAL  ::  LSECONDARY2PASTURE, LPASTURE2SECONDARY
   LOGICAL  ::  LPRIMARY,LSECONDARY,LCROP,LPASTURE
   LOGICAL  ::  INFO = .FALSE.
   LOGICAL  ::  NOPASTURE = .FALSE.

   NAMELIST /timeinfo/ year_start, year_end, ires, workingPath

   ! Read time info from namelist
   ires       = 63
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
   if (year_end < year_start) then
     write(*,*) 'Bad combination of namelist parameters'
     stop
   endif

   year_start = year_start - 1            ! To be adjusted in the time loop
   nyear      = year_end - year_start

   write(*,*) 'Work path: ',trim(workingPath)
   write(*,*) 'Processing years: ',year_start+1,year_end
   NULLIFY(MAX_EXCESS)
   DO ITIME = 1,nyear  ! open time loop
      IF (year_start + ITIME .LE. 999) THEN
         WRITE(INPUT_YEAR,'(A,i3)') '0' , year_start + ITIME
      ELSE
         WRITE(INPUT_YEAR,'(i4)') year_start + ITIME
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! INPUT-FILE: LAND COVER STATES !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FNAM = trim(workingPath) // trim(INPUT_FILES1) // 'LUH_states_' // TRIM(INPUT_YEAR) // '.nc'

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
          ELSE IF (NATTTYP(I,II).EQ.4) THEN
            IRET = nf_get_att_int(NCID,I,FATTNAM(I,II),NATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,NATT',I,TRIM(FATTNAM(I,II)),NATT(I,II)
          ELSE IF ((NATTTYP(I,II).EQ.5) .OR. (NATTTYP(I,II).EQ.6)) THEN
            IRET = nf_get_att_double(NCID,I,FATTNAM(I,II),RATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,RATT',I,TRIM(FATTNAM(I,II)),RATT(I,II)
          ELSE
            WRITE (*,*) 'Type of attribute not prepared in the code, failed to read grid',I,II
            STOP
          END IF
        END DO
      END DO
! get data
      ISTATES = 0
      ALLOCATE(FSTATENAMES(NVAR))
      FSTATENAMES(:) = ''
      DO I = 1,NVAR
         IF (.NOT. FVARNAM(I) == 'lon' .AND. .NOT. FVARNAM(I) == 'lat' .AND. .NOT. FVARNAM(I) == 'time') THEN
            IF (.NOT. ALLOCATED(RDATA)) THEN
               ALLOCATE(RDATA(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),1))
               ALLOCATE(RSTATES(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),NVAR))
            END IF
            ISTATES = ISTATES + 1
            FSTATENAMES(ISTATES) = FVARNAM(I)
            IRET = nf_get_var_double(NCID,I,RDATA)
            CALL check_err(IRET)
            RSTATES(:,:,I) = RDATA(:,:,1)
         END IF
      END DO
      DEALLOCATE(RDATA)
      IF (ISTATES /= NSTATES) STOP 'Number of states wrong!'
! close the input-file states
      IRET = nf_close(NCID)
      CALL check_err(IRET)

! check number of longitudes and latitudes
      LLON = .false.
      LLAT = .false.
      DO I = 1,NDIM
         IF (FDIMNAM(I) == 'lon') THEN
            IF (NDIMLEN(I) == NLON) LLON = .true.
         ELSEIF (FDIMNAM(I) == 'lat') THEN
            IF (NDIMLEN(I) == NLAT) LLAT = .true.
         END IF
      END DO

      IF (.NOT. LLON) STOP 'Dimension lon not found or number of longitudes not correct'
      IF (.NOT. LLAT) STOP 'Dimension lat not found or number of latitudes not correct'

      LPRIMARY = .false.
      LSECONDARY = .false.
      LCROP = .false.
      LPASTURE = .false.

      DO I = 1,NVAR
         ! 'gothr' is primary vegetation
         IF (FVARNAM(I) == 'gothr') THEN
            LPRIMARY = .true.
            IPRIMARY = I
         ! 'gsecd' is secondary vegetation
         ELSE IF (FVARNAM(I) == 'gsecd') THEN
            LSECONDARY = .true.
            ISECONDARY = I
         ELSE IF (FVARNAM(I) == 'gcrop') THEN
            LCROP = .true.
            ICROP = I
         ELSE IF (FVARNAM(I) == 'gpast') THEN
            LPASTURE = .true.
            IPASTURE = I
         END IF
      END DO

      IF (.NOT. LPRIMARY) STOP 'Variable gothr (primary) not found'
      IF (.NOT. LSECONDARY) STOP 'Variable gsecd (secondary) not found'
      IF (.NOT. LCROP) STOP 'Variable gcrop not found'
      IF (.NOT. LPASTURE) STOP 'Variable gpast not found'

! deallocate
      DEALLOCATE(NDIMLEN,FDIMNAM,FVARNAM,NVARTYP,NVARDIM,NVARDIMID,NVARATT)
      DEALLOCATE(FATTNAM,NATTTYP,NATTDIM,FATT,NATT,RATT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! INPUT-FILE TRANSITIONS!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FNAM = trim(workingPath) // trim(INPUT_FILES2) // 'LUH_transitions_' // TRIM(INPUT_YEAR) // '.nc'

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
        IF (INFO) WRITE (*,*) 'I,FVARNAM,NVARTYP,NVARDIM,NVARDIMID,NVARATT',&
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
          ELSE IF (NATTTYP(I,II).EQ.4) THEN
            IRET = nf_get_att_int(NCID,I,FATTNAM(I,II),NATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,NATT',I,TRIM(FATTNAM(I,II)),NATT(I,II)
          ELSE IF ((NATTTYP(I,II).EQ.5) .OR. (NATTTYP(I,II).EQ.6)) THEN
            IRET = nf_get_att_double(NCID,I,FATTNAM(I,II),RATT(I,II))
            CALL check_err(IRET)
            IF (INFO) WRITE (*,*) 'I,FATTNAM,RATT',I,TRIM(FATTNAM(I,II)),RATT(I,II)
          ELSE
            WRITE (*,*) 'Type of attribute not prepared in the code, failed to read grid',I,II
            STOP
          END IF
        END DO
      END DO
! get data
      ITRANS = 0
      DO I = 1,NVAR
         IF (FVARNAM(I) == 'lon' ) THEN
            ALLOCATE(RVARLON(NDIMLEN(NVARDIMID(I,1))))
            IRET = nf_get_var_double(NCID,I,RVARLON)
            CALL check_err(IRET)
         ELSE IF (FVARNAM(I) == 'lat') THEN
            ALLOCATE(RVARLAT(NDIMLEN(NVARDIMID(I,1))))
            IRET = nf_get_var_double(NCID,I,RVARLAT)
            CALL check_err(IRET)
         ELSE IF (FVARNAM(I) == 'time') THEN
            ALLOCATE(RVARTIME(NDIMLEN(NVARDIMID(I,1))))
            IRET = nf_get_var_double(NCID,I,RVARTIME)
            CALL check_err(IRET)
         ELSE
            IF (.NOT. ALLOCATED(RDATA)) THEN
               ALLOCATE(RDATA(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),1))
               ALLOCATE(RTRANS(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),NVAR))
               ALLOCATE(RTRANS_REL(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),NVAR))
               ALLOCATE(RTRANS_OUT(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),NTRANS_OUT))
               ALLOCATE(ROUTPUT(NDIMLEN(NVARDIMID(I,1)),NDIMLEN(NVARDIMID(I,2)),1))
               RTRANS(:,:,:) = 0.
               RTRANS_REL(:,:,:) = 0.
               RTRANS_OUT(:,:,:) = 0.
            END IF
            ITRANS = ITRANS + 1
            IRET = nf_get_var_double(NCID,I,RDATA)
            CALL check_err(IRET)
            RTRANS(:,:,I) = RDATA(:,:,1)
         END IF
      END DO
      DEALLOCATE(RDATA)
      IF (ITRANS /= NTRANS) STOP'Number of transitions wrong!'
! close the input-file transitions
      IRET = nf_close(NCID)
      CALL check_err(IRET)
! check number of longitudes and latitudes
      LLON = .false.
      LLAT = .false.
      DO I = 1,NDIM
         IF (FDIMNAM(I) == 'lon') THEN
            NDIMINID(1) = I
            IF (NDIMLEN(I) == NLON) LLON = .true.
         ELSE IF (FDIMNAM(I) == 'lat') THEN
            NDIMINID(2) = I
            IF (NDIMLEN(I) == NLAT) LLAT = .true.
         ELSE IF (FDIMNAM(I) == 'time') THEN
            NDIMINID(3) = I
         END IF
      END DO

      IF (.NOT. LLON) STOP 'Dimension lon not found or number of longitudes not correct'
      IF (.NOT. LLAT) STOP 'Dimension lat not found or number of latitudes not correct'

      LCROP2PASTURE = .false.
      LPASTURE2CROP = .false.
      LPRIMARY2PASTURE = .false.
      LPRIMARY2CROP = .false.
      LSECONDARY2CROP = .false.
      LCROP2SECONDARY = .false.
      LSECONDARY2PASTURE = .false.
      LPASTURE2SECONDARY = .false.

      DO I = 1,NVAR
         IF (FVARNAM(I) == 'crop2pasture') THEN
            LCROP2PASTURE = .true.
            ICROP2PASTURE = I
         ELSE IF (FVARNAM(I) == 'pasture2crop') THEN
            LPASTURE2CROP = .true.
            IPASTURE2CROP = I
         ELSE IF (FVARNAM(I) == 'primary2pasture') THEN
            LPRIMARY2PASTURE = .true.
            IPRIMARY2PASTURE = I
         ELSE IF (FVARNAM(I) == 'primary2crop') THEN
            LPRIMARY2CROP = .true.
            IPRIMARY2CROP = I
         ELSE IF (FVARNAM(I) == 'secondary2crop') THEN
            LSECONDARY2CROP = .true.
            ISECONDARY2CROP = I
         ELSE IF (FVARNAM(I) == 'crop2secondary') THEN
            LCROP2SECONDARY = .true.
            ICROP2SECONDARY = I
         ELSE IF (FVARNAM(I) == 'secondary2pasture') THEN
            LSECONDARY2PASTURE = .true.
            ISECONDARY2PASTURE = I
         ELSE IF (FVARNAM(I) == 'pasture2secondary') THEN
            LPASTURE2SECONDARY = .true.
            IPASTURE2SECONDARY = I
         END IF
      END DO
      IF (.NOT. LCROP2PASTURE) STOP 'Variable crop2pasture not found'
      IF (.NOT. LPASTURE2CROP) STOP 'Variable pasture2crop not found'
      IF (.NOT. LPRIMARY2PASTURE) STOP 'Variable primary2pasture not found'
      IF (.NOT. LPRIMARY2CROP) STOP 'Variable primary2crop not found'
      IF (.NOT. LSECONDARY2CROP) STOP 'Variable secondary2crop not found'
      IF (.NOT. LCROP2SECONDARY) STOP 'Variable crop2secondary not found'
      IF (.NOT. LSECONDARY2PASTURE) STOP 'Variable secondary2pasture not found'
      IF (.NOT. LPASTURE2SECONDARY) STOP 'Variable pasture2secondary not found'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Change from absolute transitions to relative transitions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (.NOT. ASSOCIATED(MAX_EXCESS)) THEN
         ALLOCATE(MAX_EXCESS(NVAR,4))
         MAX_EXCESS(:,:) = 0._dp
      ENDIF
      DO J = 1,NLON
         DO K = 1,NLAT
            IF (RSTATES(J,K,ICROP) > 0._dp) THEN
               RTRANS_REL(J,K,ICROP2PASTURE)   = MAX(0._dp,MIN(1._dp,RTRANS(J,K,ICROP2PASTURE)   / RSTATES(J,K,ICROP)))
               RTRANS_REL(J,K,ICROP2SECONDARY) = MAX(0._dp,MIN(1._dp,RTRANS(J,K,ICROP2SECONDARY) / RSTATES(J,K,ICROP)))
               IF (RTRANS(J,K,ICROP2PASTURE) / RSTATES(J,K,ICROP) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,ICROP2PASTURE)/RSTATES(J,K,ICROP) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'CROP 2 PASTURE     :',J,K,RTRANS(J,K,ICROP2PASTURE),RSTATES(J,K,ICROP),EXCESS
                 IF (MAX_EXCESS(ICROP2PASTURE,1) <   EXCESS) &
                     MAX_EXCESS(ICROP2PASTURE,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
               IF (RTRANS(J,K,ICROP2SECONDARY) / RSTATES(J,K,ICROP) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,ICROP2SECONDARY)/RSTATES(J,K,ICROP) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'CROP 2 SECONDARY   :',J,K,RTRANS(J,K,ICROP2SECONDARY),RSTATES(J,K,ICROP),EXCESS
                 IF (MAX_EXCESS(ICROP2SECONDARY,1) <   EXCESS) &
                     MAX_EXCESS(ICROP2SECONDARY,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
            ELSE
               RTRANS_REL(J,K,ICROP2PASTURE)   = 0._dp
               RTRANS_REL(J,K,ICROP2SECONDARY) = 0._dp
            END IF
            IF (RSTATES(J,K,IPASTURE) > 0._dp) THEN
               RTRANS_REL(J,K,IPASTURE2CROP)      = MAX(0._dp,MIN(1._dp,RTRANS(J,K,IPASTURE2CROP)      / RSTATES(J,K,IPASTURE)))
               RTRANS_REL(J,K,IPASTURE2SECONDARY) = MAX(0._dp,MIN(1._dp,RTRANS(J,K,IPASTURE2SECONDARY) / RSTATES(J,K,IPASTURE)))
               IF (RTRANS(J,K,IPASTURE2CROP) / RSTATES(J,K,IPASTURE) > 1.) THEN
                 EXCESS = RTRANS(J,K,IPASTURE2CROP)/RSTATES(J,K,IPASTURE) - 1._dp
                  IF (INFO) &
                  & WRITE (*,'(A,2i4,3g16.8)')'PASTURE 2 CROP     :',J,K,RTRANS(J,K,IPASTURE2CROP),RSTATES(J,K,IPASTURE),EXCESS
                 IF (MAX_EXCESS(IPASTURE2CROP,1) <   EXCESS) &
                     MAX_EXCESS(IPASTURE2CROP,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
               IF (RTRANS(J,K,IPASTURE2SECONDARY) / RSTATES(J,K,IPASTURE) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,IPASTURE2SECONDARY)/RSTATES(J,K,IPASTURE) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'PASTURE 2 SECONDARY:',J,K,RTRANS(J,K,IPASTURE2SECONDARY),RSTATES(J,K,IPASTURE),EXCESS
                 IF (MAX_EXCESS(IPASTURE2SECONDARY,1) <   EXCESS) &
                     MAX_EXCESS(IPASTURE2SECONDARY,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
            ELSE
               RTRANS_REL(J,K,IPASTURE2CROP)      = 0._dp
               RTRANS_REL(J,K,IPASTURE2SECONDARY) = 0._dp
            END IF
            IF (RSTATES(J,K,IPRIMARY) > 0._dp) THEN
               RTRANS_REL(J,K,IPRIMARY2PASTURE) = MAX(0._dp,MIN(1._dp,RTRANS(J,K,IPRIMARY2PASTURE) / RSTATES(J,K,IPRIMARY)))
               RTRANS_REL(J,K,IPRIMARY2CROP)    = MAX(0._dp,MIN(1._dp,RTRANS(J,K,IPRIMARY2CROP)    / RSTATES(J,K,IPRIMARY)))
               IF (RTRANS(J,K,IPRIMARY2PASTURE) / RSTATES(J,K,IPRIMARY) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,IPRIMARY2PASTURE)/RSTATES(J,K,IPRIMARY) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'PRIMARY 2 PASTURE  :',J,K,RTRANS(J,K,IPRIMARY2PASTURE),RSTATES(J,K,IPRIMARY),EXCESS
                 IF (MAX_EXCESS(IPRIMARY2PASTURE,1) <   EXCESS) &
                     MAX_EXCESS(IPRIMARY2PASTURE,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
               IF (RTRANS(J,K,IPRIMARY2CROP) / RSTATES(J,K,IPRIMARY) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,IPRIMARY2CROP)/RSTATES(J,K,IPRIMARY) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'PRIMARY 2 CROP     :',J,K,RTRANS(J,K,IPRIMARY2CROP),RSTATES(J,K,IPRIMARY),EXCESS
                 IF (MAX_EXCESS(IPRIMARY2CROP,1) <   EXCESS) &
                     MAX_EXCESS(IPRIMARY2CROP,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
            ELSE
               RTRANS_REL(J,K,IPRIMARY2PASTURE) = 0._dp
               RTRANS_REL(J,K,IPRIMARY2CROP)    = 0._dp
            END IF
            IF (RSTATES(J,K,ISECONDARY) > 0._dp) THEN
               RTRANS_REL(J,K,ISECONDARY2PASTURE) = MAX(0._dp,MIN(1._dp,RTRANS(J,K,ISECONDARY2PASTURE) / RSTATES(J,K,ISECONDARY)))
               RTRANS_REL(J,K,ISECONDARY2CROP)    = MAX(0._dp,MIN(1._dp,RTRANS(J,K,ISECONDARY2CROP)    / RSTATES(J,K,ISECONDARY)))
               IF (RTRANS(J,K,ISECONDARY2PASTURE) / RSTATES(J,K,ISECONDARY) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,ISECONDARY2PASTURE)/RSTATES(J,K,ISECONDARY) - 1._dp
                 IF (INFO) WRITE (*,'(A,2i4,3g16.8)')'SECONDARY 2 PASTURE:' &
                  & ,J,K,RTRANS(J,K,ISECONDARY2PASTURE),RSTATES(J,K,ISECONDARY),EXCESS
                 IF (MAX_EXCESS(ISECONDARY2PASTURE,1) <   EXCESS) &
                     MAX_EXCESS(ISECONDARY2PASTURE,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
               IF (RTRANS(J,K,ISECONDARY2CROP) / RSTATES(J,K,ISECONDARY) > 1._dp) THEN
                 EXCESS = RTRANS(J,K,ISECONDARY2CROP)/RSTATES(J,K,ISECONDARY) - 1._dp
                 IF (INFO) &
                 & WRITE (*,'(A,2i4,3g16.8)')'SECONDARY 2 CROP   :',J,K,RTRANS(J,K,ISECONDARY2CROP),RSTATES(J,K,ISECONDARY),EXCESS
                 IF (MAX_EXCESS(ISECONDARY2CROP,1) <   EXCESS) &
                     MAX_EXCESS(ISECONDARY2CROP,:) = (/EXCESS,REAL(ITIME+year_start,dp),REAL(J,dp),REAL(K,dp)/)
               ENDIF
            ELSE
               RTRANS_REL(J,K,ISECONDARY2PASTURE) = 0._dp
               RTRANS_REL(J,K,ISECONDARY2CROP)    = 0._dp
            END IF
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! COMPUTE NEW TRANSITION MATRIX !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(FTRANSNAMES(NTRANS_OUT))
      RTRANS_OUT(:,:,:) = 0.
      RTRANS_OUT(:,:,1) = RTRANS_REL(:,:,ICROP2PASTURE)
      FTRANSNAMES(1) = 'crop2pasture'
      RTRANS_OUT(:,:,2) = RTRANS_REL(:,:,IPASTURE2CROP)
      FTRANSNAMES(2) = 'pasture2crop'
      DO I = 1,NLON
         DO J = 1,NLAT
            IF (RTRANS_REL(I,J,IPRIMARY2CROP) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,ISECONDARY2CROP) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,IPRIMARY2CROP) <= 1. + EPSILON(1.) .AND. &
                RTRANS_REL(I,J,ISECONDARY2CROP) <= 1. + EPSILON(1.) .AND. &
                RSTATES(I,J,IPRIMARY) >= -EPSILON(1.) .AND. &
                RSTATES(I,J,ISECONDARY) >= -EPSILON(1.) .AND. &
                RSTATES(I,J,IPRIMARY) <= 1. + EPSILON(1.) .AND. &
                RSTATES(I,J,ISECONDARY) <= 1. + EPSILON(1.)) THEN
               IF ((RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY)) > EPSILON(1.)) THEN
                  RTRANS_OUT(I,J,3) = (RTRANS_REL(I,J,IPRIMARY2CROP) * RSTATES(I,J,IPRIMARY) + &
                                      RTRANS_REL(I,J,ISECONDARY2CROP) * RSTATES(I,J,ISECONDARY)) / &
                                      (RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY))
               ELSE
                  RTRANS_OUT(I,J,3) = 0.
               END IF
            ELSE
               RTRANS_OUT(I,J,3) = MISS_VALUE
            END IF
            IF (RTRANS_REL(I,J,IPRIMARY2PASTURE) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,ISECONDARY2PASTURE) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,IPRIMARY2PASTURE) <= 1. + EPSILON(1.) .AND. &
                RTRANS_REL(I,J,ISECONDARY2PASTURE) <= 1. + EPSILON(1.) .AND. &
                RSTATES(I,J,IPRIMARY) >= -EPSILON(1.) .AND. &
                RSTATES(I,J,ISECONDARY) >= -EPSILON(1.) .AND. &
                RSTATES(I,J,IPRIMARY) <= 1. + EPSILON(1.) .AND. &
                RSTATES(I,J,ISECONDARY) <= 1. + EPSILON(1.)) THEN
               IF ((RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY)) > EPSILON(1.)) THEN
                  RTRANS_OUT(I,J,4) = (RTRANS_REL(I,J,IPRIMARY2PASTURE) * RSTATES(I,J,IPRIMARY) + &
                                      RTRANS_REL(I,J,ISECONDARY2PASTURE) * RSTATES(I,J,ISECONDARY)) / &
                                      (RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY))
               ELSE
                  RTRANS_OUT(I,J,4) = RTRANS_REL(I,J,IPRIMARY2PASTURE) + RTRANS_REL(I,J,ISECONDARY2PASTURE)
               END IF
            ELSE
               RTRANS_OUT(I,J,4) = MISS_VALUE
            END IF
            IF (RTRANS_REL(I,J,ICROP2SECONDARY) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,ICROP2SECONDARY) <= 1. + EPSILON(1.)) THEN
               RTRANS_OUT(I,J,5) = RTRANS_REL(I,J,ICROP2SECONDARY)
            ELSE
               RTRANS_OUT(I,J,5) = MISS_VALUE
            END IF
            IF (RTRANS_REL(I,J,IPASTURE2SECONDARY) >= -EPSILON(1.) .AND. &
                RTRANS_REL(I,J,IPASTURE2SECONDARY) <= 1. + EPSILON(1.)) THEN
               RTRANS_OUT(I,J,6) = RTRANS_REL(I,J,IPASTURE2SECONDARY)
            ELSE
               RTRANS_OUT(I,J,6) = MISS_VALUE
            END IF
         END DO
      END DO
      FTRANSNAMES(3) = 'natural2crop'
      FTRANSNAMES(4) = 'natural2pasture'
      FTRANSNAMES(5) = 'crop2natural'
      FTRANSNAMES(6) = 'pasture2natural'
!!!!!!!!!!!!!!!!!!
!!! no pasture !!!
!!!!!!!!!!!!!!!!!!
      IF (NOPASTURE) THEN
         DO I = 1,NLON
            DO J = 1,NLAT
               RNAT_TO_NAT = 1. - RTRANS_OUT(I,J,3) - RTRANS_OUT(I,J,4)
               RSTATE_NATURAL_NEXT = RNAT_TO_NAT * (RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY)) + &
                                     RTRANS_OUT(I,J,5) * RSTATES(I,J,ICROP) + &
                                     RTRANS_OUT(I,J,6) * RSTATES(I,J,IPASTURE)
               RSTATE_MANAGED_NOW = RSTATES(I,J,ICROP) + RSTATES(I,J,IPASTURE)
               RTRANS_NATURAL_MANAGED = 1. - RNAT_TO_NAT
               IF (RSTATE_MANAGED_NOW > EPSILON(1.)) THEN
                  RTRANS_MANAGED_NATURAL = (RSTATE_NATURAL_NEXT - &
                                            RNAT_TO_NAT * (RSTATES(I,J,IPRIMARY) + RSTATES(I,J,ISECONDARY))) / &
                                            RSTATE_MANAGED_NOW
               ELSE
                  RTRANS_MANAGED_NATURAL = 0.
               END IF
               RTRANS_OUT(I,J,1) = 0.
               RTRANS_OUT(I,J,2) = 0.
               RTRANS_OUT(I,J,3) = RTRANS_NATURAL_MANAGED ! replace transition "natural to crop"
               RTRANS_OUT(I,J,4) = 0.
               RTRANS_OUT(I,J,5) = RTRANS_MANAGED_NATURAL ! replace transition "crop to natural"
               RTRANS_OUT(I,J,6) = 0.
               RSTATES(I,J,ICROP) = RSTATE_MANAGED_NOW
               RSTATES(I,J,IPASTURE) = 0.
            END DO
         END DO
      END IF ! nopasture
!!!!!!!!!!!!!!!!!!!!!!!!!
!!! OUTPUT TRANSITIONS!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (NOPASTURE) THEN
        FNAM = trim(workingPath) // trim(OUTPUT_FILES) // 'LUH_transitions_nopasture_' // TRIM(INPUT_YEAR) // '.nc'
      ELSE
        FNAM = trim(workingPath) // trim(OUTPUT_FILES) // 'LUH_transitions_' // TRIM(INPUT_YEAR) // '.nc'
      END IF

! open the output-file
      WRITE (*,*) 'Output to: ',TRIM(FNAM)
      CALL check_err(nf_create(FNAM,nf_clobber,NCIDOUT))
! define dimensions also as variables (same as the dimensions in the transitions input-file)
      NDIMLEN(:) = (/ NDIMLEN(NDIMINID(1)), NDIMLEN(NDIMINID(2)), NF_UNLIMITED /) ! dim: lon, lat, time
      DO I = 1,NDIM
         CALL check_err(nf_def_dim(NCIDOUT,FDIMNAM(NDIMINID(I)),NDIMLEN(I), NDIMOUTID(I)))
         CALL check_err(nf_def_var(NCIDOUT,FDIMNAM(NDIMINID(I)),6,1,NDIMOUTID(I),NVAROUTID(I)))
         CALL PUTATT(NVAR,NDIMINID(I),NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(I),NATT,RATT,INFO)
      END DO
! define non-dimension variables
      DO I = 1,NTRANS_OUT
        CALL check_err(nf_def_var(NCIDOUT,FTRANSNAMES(I),6,3,NDIMOUTID,NVAROUTID(I+3)))
        CALL PUTATT(NVAR,I+3,NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(I+3),NATT,RATT,INFO)
      ENDDO
      CALL check_err(nf_enddef(NCIDOUT))
! put longitudes, latitudes and time
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(1),RVARLON))
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(2),RVARLAT))
      CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(3),(/1/),(/1/),RVARTIME))
! put data
      DO I = 1,NTRANS_OUT
         ROUTPUT(:,:,1) = RTRANS_OUT(:,:,I)
         CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(I+3),(/1,1,1/),(/NDIMLEN(1),NDIMLEN(2),1/),ROUTPUT(:,:,1)))
      END DO

! close the output-file
      CALL check_err(nf_close(NCIDOUT))

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
      ALLOCATE(RSTATES_FLIP(NLON+2,NLAT+2))
      ALLOCATE(RTRANS_OUT_GAUSS(NGLON,NGLAT,NTRANS_OUT))
      ALLOCATE(ROUTPUT_GAUSS(NGLON,NGLAT,1))
      RTRANS_OUT_GAUSS(:,:,:) = 0.

      DO II = 1,NTRANS_OUT
         IF (INFO) WRITE(*,'(A,i1)')  'transition nr ', II
         RTRANS_FLIP(:,:) = 0.
         RSTATES_FLIP(:,:) = 0.
         IF (II == 1) ISTATES_OUT = ICROP
         IF (II == 2) ISTATES_OUT = IPASTURE
         IF (II == 3) ISTATES_OUT = IPRIMARY
         IF (II == 4) ISTATES_OUT = IPRIMARY
         IF (II == 5) ISTATES_OUT = ICROP
         IF (II == 6) ISTATES_OUT = IPASTURE
         DO J = 1,NLAT
            DO I = 1,NLON/2
               RTRANS_FLIP(I + 1,J + 1) = RTRANS_OUT(I + NLON/2,J,II)
               RSTATES_FLIP(I + 1,J + 1) = RSTATES(I + NLON/2,J,ISTATES_OUT)
               IF (II == 3 .OR. II == 4) THEN
                  IF (RSTATES_FLIP(I + 1,J + 1) >= -EPSILON(1.) .AND. &
                      RSTATES_FLIP(I + 1,J + 1) <= 1. + EPSILON(1.) .AND. &
                      RSTATES(I + NLON/2,J,ISECONDARY) >= -EPSILON(1.) .AND. &
                      RSTATES(I + NLON/2,J,ISECONDARY) <= 1. + EPSILON(1.)) THEN
                     RSTATES_FLIP(I + 1,J + 1) = RSTATES_FLIP(I + 1,J + 1) + RSTATES(I + NLON/2,J,ISECONDARY)
                  ELSE
                     RSTATES_FLIP(I + 1,J + 1) = MISS_VALUE
                  END IF
               END IF
            END DO
            DO I = NLON/2 + 1,NLON
               RTRANS_FLIP(I + 1,J + 1) = RTRANS_OUT(I - NLON/2,J,II)
               RSTATES_FLIP(I + 1,J + 1) = RSTATES(I - NLON/2,J,ISTATES_OUT)
               IF (II == 3 .OR. II == 4) THEN
                  IF (RSTATES_FLIP(I + 1,J + 1) >= -EPSILON(1.) .AND. &
                      RSTATES_FLIP(I + 1,J + 1) <= 1. + EPSILON(1.) .AND. &
                      RSTATES(I - NLON/2,J,ISECONDARY) >= -EPSILON(1.) .AND. &
                      RSTATES(I - NLON/2,J,ISECONDARY) <= 1. + EPSILON(1.)) THEN
                     RSTATES_FLIP(I + 1,J + 1) = RSTATES_FLIP(I + 1,J + 1) + RSTATES(I - NLON/2,J,ISECONDARY)
                  ELSE
                     RSTATES_FLIP(I + 1,J + 1) = MISS_VALUE
                  END IF
               END IF
            END DO
            RTRANS_FLIP(1,J + 1) = RTRANS_FLIP(NLON + 1,J + 1)
            RTRANS_FLIP(NLON + 2, J + 1) = RTRANS_FLIP(2,J + 1)
            RSTATES_FLIP(1,J + 1) = RSTATES_FLIP(NLON + 1,J + 1)
            RSTATES_FLIP(NLON + 2, J + 1) = RSTATES_FLIP(2,J + 1)
         END DO
         DO I = 1,NLON + 2
            RTRANS_FLIP(I,1) = RTRANS_FLIP(I,2)
            RTRANS_FLIP(I,NLAT + 2) = RTRANS_FLIP(I,NLAT + 1)
            RSTATES_FLIP(I,1) = RSTATES_FLIP(I,2)
            RSTATES_FLIP(I,NLAT + 2) = RSTATES_FLIP(I,NLAT + 1)
         END DO
         IF (INFO) WRITE (*,*) 'call intpol2'
         CALL INTPOL2(INFO,NLON+2,NLAT+2,RTRANS_FLIP,RSTATES_FLIP,XILON,XILAT,NGLON,NGLAT,XGLON,XGLAT,MISS_VALUE, &
                      RTRANS_OUT_GAUSS(:,:,II))
!       RTRANS_OUT_GAUSS(:,:,II) = 0.0
      END DO
      IF (MAXVAL(RTRANS_OUT_GAUSS(:,:,:)) > 1.0000001) STOP 'transition larger than 1'
      IF (MINVAL(RTRANS_OUT_GAUSS(:,:,:)) < -0.0000001) STOP 'transition smaller than 1'
      RTRANS_OUT_GAUSS(:,:,:) = MIN(1._dp-EPSILON(1._dp),MAX(0._dp,RTRANS_OUT_GAUSS(:,:,:)))
      IF (NOPASTURE) THEN
         FNAM = trim(workingPath) // trim(OUTPUT_FILES) // &
           & 'LUH_transitions_nopasture_T' // TRIM(FRES) // '_' // TRIM(INPUT_YEAR) // '.nc'
      ELSE
         FNAM = trim(workingPath) // trim(OUTPUT_FILES) // 'LUH_transitions_T' // TRIM(FRES) // '_' // TRIM(INPUT_YEAR) // '.nc'
      ENDIF
! switch latitudes
      DO J = 1,NGLAT
         XGLAT_SWITCH(J) = XGLAT(NGLAT + 1 -J)
      END DO
! open the output-file
      WRITE (*,*) 'Output to: ',TRIM(FNAM)
      CALL check_err(nf_create(FNAM,nf_clobber,NCIDOUT))
! define dimensions as dimensions and variables (same as the dimensions in the transitions input-file)
      NDIMLEN(1:2) = (/NGLON,NGLAT/) ! lon, lat, time (not changed)
      DO I = 1,NDIM
         CALL check_err(nf_def_dim(NCIDOUT,FDIMNAM(NDIMINID(I)),NDIMLEN(I),NDIMOUTID(I)))
         CALL check_err(nf_def_var(NCIDOUT,FVARNAM(NDIMINID(I)),6,1,NDIMOUTID(I:I),NVAROUTID(I)))
         CALL PUTATT(NVAR,NDIMINID(I),NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(I),NATT,RATT,INFO)
      END DO
! define other variables
      DO I = 1,NTRANS_OUT
         CALL check_err(nf_def_var(NCIDOUT,FTRANSNAMES(I),6,3,NDIMOUTID,NVAROUTID(I+3)))
         CALL PUTATT(NVAR,I+3,NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVAROUTID(I+3),NATT,RATT,INFO)
      ENDDO
      CALL check_err(nf_enddef(NCIDOUT))
! put longitudes, latitudes and time
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(1),XGLON))
      CALL check_err(nf_put_var_double (NCIDOUT,NVAROUTID(2),XGLAT_SWITCH))
      CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(3),(/1/),(/1/),RVARTIME))
! put data

      DO I = 1,NTRANS_OUT
         DO J = 1,NGLAT
            ROUTPUT_GAUSS(:,J,1) = RTRANS_OUT_GAUSS(:,NGLAT + 1 - J,I)
         END DO
         CALL check_err(nf_put_vara_double(NCIDOUT,NVAROUTID(I+3),(/1,1,1/),(/NGLON,NGLAT,1/),ROUTPUT_GAUSS(:,:,1)))
      END DO

! close the output-file
      CALL check_err(nf_close(NCIDOUT))

! deallocate arrays
      DEALLOCATE(NDIMLEN,FDIMNAM,FVARNAM,NVARTYP,NVARDIM,NVARDIMID,NVARATT)
      DEALLOCATE(FATTNAM,NATTTYP,NATTDIM,FATT,NATT,RATT,RVARLON,RVARLAT,RVARTIME)
      DEALLOCATE(FTRANSNAMES,FSTATENAMES)
      DEALLOCATE(RSTATES,RTRANS,RTRANS_REL)
      DEALLOCATE(RTRANS_OUT,ROUTPUT)
      DEALLOCATE(XILAT,XILON,XGLAT,XGLON,XGLAT_SWITCH)
      DEALLOCATE(RTRANS_FLIP,RSTATES_FLIP,RTRANS_OUT_GAUSS,ROUTPUT_GAUSS)

   END DO  ! close time loop

! print error "statistics"
   WRITE(*,'(A,g16.8,3i5)') 'Max excess PASTURE   2 CROP     :',MAX_EXCESS(IPASTURE2CROP,1)     , &
                                                            INT(MAX_EXCESS(IPASTURE2CROP,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess PASTURE   2 SECONDARY:',MAX_EXCESS(IPASTURE2SECONDARY,1), &
                                                            INT(MAX_EXCESS(IPASTURE2SECONDARY,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess CROP      2 PASTURE  :',MAX_EXCESS(ICROP2PASTURE,1)     , &
                                                            INT(MAX_EXCESS(ICROP2PASTURE,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess CROP      2 SECONDARY:',MAX_EXCESS(ICROP2SECONDARY,1)   , &
                                                            INT(MAX_EXCESS(ICROP2SECONDARY,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess PRIMARY   2 CROP     :',MAX_EXCESS(IPRIMARY2CROP,1)     , &
                                                            INT(MAX_EXCESS(IPRIMARY2CROP,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess PRIMARY   2 PASTURE  :',MAX_EXCESS(IPRIMARY2PASTURE,1)  , &
                                                            INT(MAX_EXCESS(IPRIMARY2PASTURE,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess SECONDARY 2 PASTURE  :',MAX_EXCESS(ISECONDARY2PASTURE,1), &
                                                            INT(MAX_EXCESS(ISECONDARY2PASTURE,2:4))
   WRITE(*,'(A,g16.8,3i5)') 'Max excess SECONDARY 2 CROP     :',MAX_EXCESS(ISECONDARY2CROP,1)   , &
                                                            INT(MAX_EXCESS(ISECONDARY2CROP,2:4))
   DEALLOCATE(MAX_EXCESS)

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

   SUBROUTINE INTPOL2(INFO,NXM,NYM,FIELDM,WEIGHT,XM,YM,NX,NY,XX,YY,FILL_VALUE,FIELD)
   ! Flux conserving transformation of a global field with a high resolution to a lower resolution grid.
   ! Weighting by fractional cover.
   ! The incoming field is assumed to range from the South Pole to the North Pole.
   ! This routine is designed for pre- and post-processing purposes It is not optimized in terms of computational efficiency.
   ! From the tested version /pf/m/m212070/save/C4MIP/OCMIP/lidsst8.f90
   ! ATTENTION: If the outgoing grid is finer than the incoming grid the result will be erroneous without warning
   ! author: Thomas Raddatz 2003 (modified for weighting November 2009)
   IMPLICIT NONE
   LOGICAL,INTENT(IN)   ::  INFO
   INTEGER,INTENT(IN)   ::  NXM,NYM
   INTEGER,INTENT(IN)   ::  NX,NY
   REAL(dp),INTENT(IN)   ::  FIELDM(NXM,NYM)
   REAL(dp),INTENT(IN)   ::  WEIGHT(NXM,NYM)
   REAL(dp),INTENT(IN)   ::  XM(NXM),YM(NYM)
   REAL(dp),INTENT(IN)   ::  XX(NX),YY(NY)
   REAL(dp),INTENT(IN)   ::  FILL_VALUE
   REAL(dp),INTENT(OUT)  ::  FIELD(NX,NY)
! local variables
   REAL(dp),PARAMETER  ::  deg_to_rad=3.1415926536/180.
   REAL(dp) :: AREA(NX,NY),XDEL,AIN,AIN1,AIN2,AIN3,AIN4,PI
   REAL(dp) :: X(NX+3),Y(NY+2)
   INTEGER :: IFLAG,IX,IY,ICASE
   LOGICAL :: speadUp = .TRUE.
! indices
   INTEGER :: I,J,II,JJ
   PI = 3.14159265359


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
!     IF (INFO) WRITE (*,*) 'longitude counter: ',II
!     IF (INFO) WRITE (*,*) II
     DO 98 JJ=2,NYM-1
       IF (FIELDM(II,JJ) > -EPSILON(1.) .AND. FIELDM(II,JJ) <= 1. + EPSILON(1.) .AND. &
           WEIGHT(II,JJ) > -EPSILON(1.) .AND. WEIGHT(II,JJ) <= 1. + EPSILON(1.)) THEN
       IFLAG=0
       AIN = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((YM(JJ) + YM(JJ+1)) *deg_to_rad/ 2.)) &
             * XDEL * 180. / PI
! assign lower left corner of the IN GRID BOX to one OUT GRID BOX
       DO I=2,NX+2
         DO J=2,NY+1
           IF (XM(II) + XM(II-1) >= X(I) + X(I-1) .AND. XM(II) + XM(II-1) < X(I) + X(I+1) .AND. &
               YM(JJ) + YM(JJ-1) >= Y(J) + Y(J-1) .AND. YM(JJ) + YM(JJ-1) < Y(J) + Y(J+1)) THEN
             IFLAG = IFLAG + 1
             IX = I
             IY = J
             IF(speadUp) exit
           END IF
         END DO
         IF (IFLAG == 1 .AND. speadUp) exit
       END DO
       IF (IFLAG == 0) STOP 'FAILED TO FIND IX,IY'
       IF (IFLAG >= 2) STOP 'FOUND MULTIPLE IX,IY'
! 4 cases:
       IFLAG = 0
       IF (XM(II) + XM(II+1) < X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) < Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX is completely in ONE OUT-GRID BOX
         IFLAG = IFLAG + 1
         ICASE = 1
       ELSE IF (XM(II) + XM(II+1) < X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) >= Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX stretches over 2 OUT-GRID BOXES in meridional direction
         IFLAG = IFLAG + 1
         ICASE = 2
       ELSE IF (XM(II) + XM(II+1) >= X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) < Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX stretches over 2 OUT-GRID BOXES in longitudinal direction
         IFLAG = IFLAG+1
         ICASE = 3
       ELSE IF (XM(II) + XM(II+1) >= X(IX) + X(IX+1) .AND. YM(JJ) + YM(JJ+1) >= Y(IY) + Y(IY+1)) THEN
! IN-GRID BOX partially covers 4 OUT-GRID BOXES
         IFLAG = IFLAG + 1
         ICASE = 4
       END IF
       IF (IX == NX+2) IX = 2
       IF (IFLAG == 0) STOP 'FAILED TO FIND ICASE'
       IF (IFLAG.GE.2) STOP 'FOUND MULTIPLE ICASE'
! add flux multiplied by the associated area to OUT-GRID
       IF (ICASE == 1) THEN
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + (WEIGHT(II,JJ) * AIN)
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN)
       ELSE IF (ICASE == 2) THEN
         IF (IY > NY) WRITE(*,*) 'IY,NY: ATTENTION (ICASE=2)!',IY,NY
         AIN1 = ABS(SIN((YM(JJ-1) + YM(JJ)) *deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) &
                * XDEL * 180. / PI
         AIN2 = ABS(SIN((YM(JJ) + YM(JJ+1)) *deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) &
                * XDEL * 180. / PI
         IF (ABS(AIN1 + AIN2 - AIN) / AIN > 1.E-04) STOP 'AREAS 2 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + (WEIGHT(II,JJ) * AIN1)
         AREA(IX-1,IY) = AREA(IX-1,IY) + (WEIGHT(II,JJ) * AIN2)
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN1)
         FIELD(IX-1,IY) = FIELD(IX-1,IY) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN2)
       ELSE IF (ICASE == 3) THEN
         AIN1 = AIN * 0.5 * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / XDEL
         AIN2 = AIN * 0.5 * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / XDEL
         IF (ABS(AIN1 + AIN2 - AIN) / AIN > 1.E-04) STOP 'AREAS 3 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + (WEIGHT(II,JJ) * AIN1)
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN1)
         IF (IX == NX+1) THEN
           AREA(1,IY-1) = AREA(1,IY-1) + (WEIGHT(II,JJ) * AIN2)
           FIELD(1,IY-1) = FIELD(1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN2)
         ELSE
           AREA(IX,IY-1) = AREA(IX,IY-1) + (WEIGHT(II,JJ) * AIN2)
           FIELD(IX,IY-1) = FIELD(IX,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN2)
         END IF
       ELSE IF (ICASE.EQ.4) THEN
         IF (IY > NY) WRITE (*,*) 'IY,NY: ACHTUNG (ICASE=4)!',IY,NY
         AIN1 = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * 90. * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * PI)
         AIN2 = ABS(SIN((YM(JJ-1) + YM(JJ)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * 90. * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * PI)
         AIN3 = ABS(SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * 90. * ABS(XM(II) + XM(II-1) - X(IX) - X(IX+1)) / (XDEL * PI)
         AIN4 = ABS(SIN((YM(JJ) + YM(JJ+1)) * deg_to_rad / 2.) - SIN((Y(IY) + Y(IY+1)) * deg_to_rad / 2.)) * &
                XDEL * 90. * ABS(XM(II) + XM(II+1) - X(IX) - X(IX+1)) / (XDEL * PI)
         IF (ABS(AIN1 + AIN2 + AIN3 + AIN4 - AIN) / AIN > 1.E-04) STOP 'AREAS 4 WRONG'
         AREA(IX-1,IY-1) = AREA(IX-1,IY-1) + (WEIGHT(II,JJ) * AIN1)
         AREA(IX-1,IY) = AREA(IX-1,IY) + (WEIGHT(II,JJ) * AIN3)
         FIELD(IX-1,IY-1) = FIELD(IX-1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN1)
         FIELD(IX-1,IY) = FIELD(IX-1,IY) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN3)
         IF (IX == NX+1) THEN
           AREA(1,IY-1) = AREA(1,IY-1) + (WEIGHT(II,JJ) * AIN2)
           AREA(1,IY) = AREA(1,IY) + (WEIGHT(II,JJ) * AIN4)
           FIELD(1,IY-1) = FIELD(1,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN2)
           FIELD(1,IY) = FIELD(1,IY) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN4)
         ELSE
           AREA(IX,IY-1) = AREA(IX,IY-1) + (WEIGHT(II,JJ) * AIN2)
           AREA(IX,IY) = AREA(IX,IY) + (WEIGHT(II,JJ) * AIN4)
           FIELD(IX,IY-1) = FIELD(IX,IY-1) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN2)
           FIELD(IX,IY) = FIELD(IX,IY) + (FIELDM(II,JJ) * WEIGHT(II,JJ) * AIN4)
         END IF
       END IF
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

END PROGRAM Calculate_relative_transitions_remap_and_flipLats
