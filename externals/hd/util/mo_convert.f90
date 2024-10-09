! mo_convert.f90 - Utilities for the conversion of HD discharges into inflows on an ocean grid.
! 
! Copyright (C) 2022, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_convert

  !
  ! Authors:
  !
  ! S. Hagemann, Hereon - Institue of Coastal Systems, January 2022, original source
  !            
  ! *** Separated from mo_interpol 

  use netcdf
  use mo_grid,         ONLY: model

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: convert_inflow, open_coupling_file, read_coupling_file


  !-----------------------------------------------------------------------------
CONTAINS

! *********************************************************************
  SUBROUTINE convert_inflow(dninp, dnout, ctitle, iconv, ybeg, yend, isep_in, isep_out, &
             mask_src, mask_target, xmiss, &
             mask_src_mapped, ix_target, iy_target, mask_on_ocean)
! *********************************************************************
!
!    Routine that uses the mapping information to transfer the inflow (e.g. discharge)
!    from an input file onto the ocean grid and the associated mouth points.
!
!    ******** Version 2.0 - June 2020
!    Programmierung und Entwicklung: Stefan Hagemann 
!    Genrealisation of program convert_hdtoocean.f90 included in script convert_hdtocean.com
!      
!    Test output is possible for a selected coordinate illog, iblog, e.g.:
!              Elbe: Lon = 8.5, Lat=54.5   Cat=12
!                    HD:   378      72
!           HD_NEMO 3.3:   233     174   after nearest neighbor interpolation
!        NEMO 3.3-coast:   230     174
!           HD_NEMO 3.6:   517     426   after nearest neighbor interpolation (EHYPE)
!
!    ******** Version 2.1 - Nov. 2021
!    Sums of inflows on in- and output grid were multiplied by the mapping masks that comprise values larger 1.
!    This bug is corrected now. 
!

    use mo_flow_inout,   ONLY: open_inflow, open_outflow, nvar, close_inflow, &
                               read_inflow, write_inflow_on_ocean, close_outflow, &
                               read_ute_inflow, set_hereon_attributes
    use mo_time,         ONLY: calc_timestep, next_timestep
    use mo_interpol,     ONLY: nsep, flow_separation
	  
    CHARACTER (LEN=*), INTENT(inout)  :: dninp   ! Input file name, general form
    CHARACTER (LEN=*), INTENT(inout)  :: dnout   ! Output file name
    CHARACTER (LEN=*), INTENT(in)     :: ctitle  ! Description of converted data
    INTEGER, INTENT(in) :: iconv    ! Method of conversion: 0 = no, 1=Daily data (def.),
                                    ! 2=monthly clim., 3=5 bgc sequence
                                    ! 4=General sequence of bgc variables and timesteps
    INTEGER, INTENT(in) :: ybeg     ! Start year
    INTEGER, INTENT(in) :: yend     ! End year
    INTEGER, INTENT(in) :: isep_in  ! Input data in separate files (No/Yes=0/1)
    INTEGER, INTENT(in) :: isep_out ! Output data in separate files (No/Yes=0/1)
!
    TYPE(model), INTENT(in)        :: mask_src    ! Source info & coordinates
    TYPE(model), INTENT(in)        :: mask_target ! Target info, e.g.
    DOUBLE PRECISION, INTENT(in)   :: xmiss    ! Missing Value, e.g. -9999.
!
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: mask_src_mapped  ! Mask of source mouths mapped to target mouth
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: ix_target        ! x-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_src%nlon,mask_src%nlat), &
        INTENT(in) :: iy_target        ! y-Indices of nearest target mouth points
    INTEGER, DIMENSION(mask_target%nlon,mask_target%nlat), &
        INTENT(in) :: mask_on_ocean    ! Mask on ocean target receiving mapping
!
!   ******* Local Variables
    INTEGER, DIMENSION(:,:), ALLOCATABLE   ::  mask_in    ! Mask of used source mouth points
    INTEGER, DIMENSION(:,:), ALLOCATABLE   ::  mask_out   ! Mask of used ocean target points
!
!   *** NETCDF variables
    INTEGER :: ierr = 0
    CHARACTER (LEN=120) :: clong 
!
    DOUBLE PRECISION, PARAMETER :: zeps = 1.E-10 
    INTEGER, PARAMETER :: ique = 0    ! Print out some debug statements for ique=1 at illog,iblog
    INTEGER, PARAMETER :: illog =378  ! Log Output for selected mouth point with lon index=illog
    INTEGER, PARAMETER :: iblog =72   ! Log Output for selected mouth point with lat index=iblog

    INTEGER :: iday, imon, iyear, iy, i, ipos_in, ipos_out
    INTEGER :: nstep_in, nstep_out
    DOUBLE PRECISION  :: yyyymmdd = 0.  ! date + time in unit [day]
    LOGICAL           :: lendofyear
!   Inflow arrays, e.g. HD
    REAL, DIMENSION(:,:), ALLOCATABLE :: friv_src
!   Target ocean grid arrays
    REAL, DIMENSION(:,:), ALLOCATABLE :: friv_on_ocean
!  
!   *** Array definition source grid, e.g. HD discharge
    ALLOCATE(friv_src(mask_src%nlon,mask_src%nlat))
    ALLOCATE(mask_in(mask_src%nlon,mask_src%nlat))
    WHERE (mask_src_mapped.GT.0.5) ; mask_in(:,:) = 1
    ELSE WHERE ; mask_in(:,:) = 0 ; END WHERE
    WRITE(*,*) "Number of Inflow points", SUM(mask_in)

!   *** Array allocation target ocean grid
    ALLOCATE(friv_on_ocean(mask_target%nlon,mask_target%nlat))
    ALLOCATE(mask_out(mask_target%nlon,mask_target%nlat))
    WHERE (mask_on_ocean.GT.0.5) ; mask_out(:,:) = 1
    ELSE WHERE ; mask_out(:,:) = 0 ; END WHERE
    WRITE(*,*) "Number of Ocean model mouth points", SUM(mask_out)
!
!   *** if isep_in/isep_out = 1, it is expected that dninp/dnout contains the string YYYY
    IF (isep_in.EQ.1) ipos_in = scan(dninp, 'YYYY')
    IF (isep_out.EQ.1) ipos_out = scan(dnout, 'YYYY')
!
!   *** Loop over years
    DO iy = ybeg, yend
      iyear = iy
      imon=1 ; iday = 0
      lendofyear = .FALSE. ; ierr = 0
      IF (isep_in.EQ.1) WRITE(dninp(ipos_in:ipos_in+3), '(I4)') iyear
      IF (isep_out.EQ.1) WRITE(dnout(ipos_out:ipos_out+3), '(I4)') iyear
!
!     *** Open inflow data and check time variable
      IF (isep_in.EQ.1 .OR. iy.EQ.ybeg) THEN
        CALL open_inflow(dninp)
        WRITE(*,*) 'Number of inflow variables: ', nvar
        nstep_in=0
        IF (mask_src%nlat.EQ.1 .AND. mask_src%nlon.EQ.376 .AND. ybeg.NE.1940) THEN    ! Utes data
          CALL calc_timestep(1, 1, imon, iyear, 1940, nstep_in)
          nstep_in = nstep_in - 1
        ENDIF
      ENDIF
!
!     *** Open discharge on ocean grid output file and write coordinates
      IF (isep_out.EQ.1 .OR. iy.EQ.ybeg) THEN
        CALL open_outflow(dnout, mask_target, iconv)
        CALL set_hereon_attributes(ctitle)
        nstep_out=0 
      ENDIF
      WRITE(*,*) 'Year ', iyear, ' nstep_in',  nstep_in, ' nstep_out',  nstep_out 
!
!     *** Read, Convert & Write within one year
      IF (iconv.LE.2) THEN
        DO WHILE (ierr.EQ.0) 
          nstep_in = nstep_in + 1
          nstep_out = nstep_out + 1
          DO i=1, nvar 
            IF (mask_src%nlat.NE.1 .OR. mask_src%nlon.NE.376) THEN
              CALL read_inflow(mask_src%nlon, mask_src%nlat, nstep_in, i, friv_src, yyyymmdd, ierr) 
              IF (nsep.GT.0) &
                 CALL flow_separation(mask_src%nlon, mask_src%nlat, friv_src)
            ELSE   ! Read Ute's data
              CALL read_ute_inflow(mask_src%nlon, nstep_in, i, xmiss, friv_src(:,1), yyyymmdd, ierr) 
            ENDIF
            IF (nstep_out.LE.2) WRITE(*,'(I3, A,F11.2,A,I4,A,I2,A,G16.6)') &
                       nstep_out, '. Input at ', yyyymmdd,' ierr=', ierr, &
                       ' Max of Variable ', i, ': ', MAXVAL(friv_src)
            IF (ierr.EQ.0) THEN
              IF(i.EQ.1) THEN
                IF (iconv.EQ.1) THEN
                  CALL next_timestep(1, iday,imon,iyear, lendofyear)   ! Next day
                  IF (iyear.NE.iy) THEN
                    WRITE(*,*) '!!! Year mismatch between loop and nextstep calculation ', iy, iyear
                    STOP 'ABBRUCH'
                  ENDIF
                ELSE IF (iconv.EQ.2) THEN
                  iday=15 
                  CALL next_timestep(2, iday,imon,iyear, lendofyear)   ! Next month
                ENDIF
                yyyymmdd = DBLE(iyear)*10000. + DBLE(imon)*100. + DBLE(iday)
!!              yyyymmdd = yyyymmdd + (hr*3600+min*60+sec)/86400.
!
              ENDIF   ! End if Variable = 1
              WHERE (mask_src_mapped(:,:).LT.0.5)
                friv_src(:,:) = 0.
              END WHERE
              IF (nstep_in.LE.5.AND.i.EQ.1) WRITE(*,*) yyyymmdd, ' -->',iyear,imon, iday,  MAXVAL(friv_src)
!    
!             ******** Transfer HD model river discharge to NEMO grid
              CALL dis_to_ocean(mask_src%nlon, mask_src%nlat, friv_src, &
                   ix_target, iy_target, mask_target%nlon, mask_target%nlat, friv_on_ocean)
              IF (nstep_in.LE.5) THEN
                WRITE(*,'(I4,A,I2,2(A,F11.2))') nstep_in, ' Var.', i, ". Sum Inflow: ", &
                     SUM(friv_src(:,:) * mask_in(:,:), ABS(friv_src-xmiss).GT.zeps), &
                     "  on Ocean: ", SUM(friv_on_ocean(:,:) * mask_out(:,:), &
                     ABS(friv_on_ocean-xmiss).GT.zeps)
              ENDIF
!
              IF (ique.GT.0) THEN
                WRITE(*,*) 'Log output discharge HD:      ', friv_src(illog, iblog), &
                    ' at ', illog, iblog 
                WRITE(*,*) 'Log output discharge HD-NEMO: ', &
                    friv_on_ocean(ix_target(illog, iblog),iy_target(illog, iblog)), &
                    ' at ', ix_target(illog, iblog),iy_target(illog, iblog)
              ENDIF
!
!             ******** Write river discharge on target ocean grid
              CALL write_inflow_on_ocean(mask_target%nlon, mask_target%nlat, &
                 nstep_out, i, friv_on_ocean, yyyymmdd)
            ENDIF
          ENDDO  ! End of Loop over Variables
          IF (lendofyear) ierr=365
        ENDDO    ! End of while loop

      ELSE IF (iconv.EQ.3) THEN  ! without time, Is this now different from iconv = 4?
        DO i=1, nvar 
          CALL read_inflow(mask_src%nlon, mask_src%nlat, 1, i, friv_src, yyyymmdd, ierr) 

          IF (ierr.EQ.0) THEN
            WHERE (mask_src_mapped(:,:).LT.0.5)
              friv_src(:,:) = 0.
            END WHERE
            WRITE(*,*) 'Variable ', i, ' -->', MINVAL(friv_src), ' - ', MAXVAL(friv_src)
!    
!           ******** Transfer HD model river discharge to NEMO grid
            CALL dis_to_ocean(mask_src%nlon, mask_src%nlat, friv_src, &
                 ix_target, iy_target, mask_target%nlon, mask_target%nlat, friv_on_ocean)
            WRITE(*,*) i, ". Sums HD: ", SUM(friv_src(:,:) * mask_in(:,:)), &
                       "  on Ocean: ", SUM(friv_on_ocean(:,:) * mask_out(:,:))
!
!           ******** Write river discharge on NEMO grid
            CALL write_inflow_on_ocean(mask_target%nlon, mask_target%nlat, &
                 1, i, friv_on_ocean, yyyymmdd)
          ENDIF
        ENDDO    ! End of loop over variables
      ELSE       ! with time
        nstep_in=0  
        ierr = 0 
        DO WHILE (ierr.EQ.0) 
          nstep_in = nstep_in + 1
          DO i=1, nvar 
            CALL read_inflow(mask_src%nlon, mask_src%nlat, nstep_in, i, friv_src, yyyymmdd, ierr) 
            IF (ierr.NE.0) THEN
              WRITE(*,*) 'End of time series (or read error) in step ', &
                          nstep_in, ' for Var: ', i, ierr
              EXIT
            ENDIF

            IF (nstep_in.LE.5 .AND. i.EQ.1) WRITE(*,*) 'Step: ', nstep_in, ierr, ' Time: ', yyyymmdd
            WHERE (mask_src_mapped(:,:).LT.0.5)
              friv_src(:,:) = 0.
            END WHERE
            IF (nstep_in.LE.5) WRITE(*,*) nstep_in, ' Variable ', i,' -->', &
                   MINVAL(friv_src), ' - ', MAXVAL(friv_src)
!    
!           ******** Transfer HD model river discharge to NEMO grid
            CALL dis_to_ocean(mask_src%nlon, mask_src%nlat, friv_src, &
                 ix_target, iy_target, mask_target%nlon, mask_target%nlat, friv_on_ocean)
            WRITE(*,*) nstep_in, ". Sums HD: ", SUM(friv_src(:,:) * mask_in(:,:)), &
                     "  on Ocean: ", SUM(friv_on_ocean(:,:) * mask_out(:,:))
!
!           ******** Write bgc flows on ocean grid
            CALL write_inflow_on_ocean(mask_target%nlon, mask_target%nlat, nstep_in, i, &
                                       friv_on_ocean, yyyymmdd)
          ENDDO    ! End of loop over variables
        ENDDO    ! End of loop over time steps within a year
      ENDIF    ! End of iconv if
      
      IF (isep_in.EQ.1) CALL close_inflow
      IF (isep_out.EQ.1) CALL close_outflow

      IF (iconv.LE.2) THEN
        WRITE(*,'(A,/,A)') "River-Discharge/BGC Inflow into ocean was written: ", TRIM(dnout)
      ELSE
        WRITE(*,'(A,/,A)') "BGC Inflow into ocean was written: ", TRIM(dnout)
      ENDIF
    ENDDO      ! End of loop over years
     
    IF (isep_in.EQ.0) CALL close_inflow
    IF (isep_out.EQ.0) CALL close_outflow
    DEALLOCATE(friv_src)
    DEALLOCATE(friv_on_ocean)

  END SUBROUTINE convert_inflow
!
! *********************************************************************
  SUBROUTINE dis_to_ocean(NL, NB, friv_src, ix_target, iy_target, nx, ny, friv_on_ocean)
! *********************************************************************
!
!   Transfer HD discharge to NEMO grid

    INTEGER, INTENT(in) :: NL, NB, nx, ny
    REAL, DIMENSION(NL,NB), INTENT(in) :: friv_src
    INTEGER, DIMENSION(NL,NB), INTENT(in) :: ix_target
    INTEGER, DIMENSION(NL,NB), INTENT(in) :: iy_target
    REAL, DIMENSION(nx,ny), INTENT(out) :: friv_on_ocean

    INTEGER :: JL, JB
!
    friv_on_ocean(:,:) = 0.
    DO JB=1, NB
    DO JL=1, NL 
    IF (ix_target(JL,JB).GT.0.5) THEN
      friv_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) = friv_on_ocean(ix_target(JL,JB), iy_target(JL,JB)) + friv_src(JL,JB)
    ENDIF
    ENDDO	  
    ENDDO	  
!
  END SUBROUTINE dis_to_ocean
!
! *********************************************************************
  SUBROUTINE open_coupling_file(dn_couple, mask_src, mask_target)
! *********************************************************************

    CHARACTER (LEN=*), INTENT(in)   :: dn_couple   ! Coupling file name
    TYPE(model), INTENT(out)        :: mask_src    ! Source grid info 
    TYPE(model), INTENT(out)        :: mask_target ! Target grid info

    INTEGER :: ncid, dimid, ierr  
!
!   ******** READ array dimensions
    ierr = nf90_open(dn_couple, NF90_NOWRITE, ncid)

    ierr = nf90_inq_dimid(ncid,'lon',dimid)
    ierr = nf90_inquire_dimension(ncid,dimid,len=mask_src%nlon)
    ierr = nf90_inq_dimid(ncid,'lat',dimid)
    ierr = nf90_inquire_dimension(ncid,dimid,len=mask_src%nlat)
    WRITE(*,*) "Coupling dimensions HD model grid: NL =" , mask_src%nlon, "  NB = ", mask_src%nlat

    ierr = nf90_inq_dimid(ncid,'x',dimid)
    ierr = nf90_inquire_dimension(ncid,dimid,len=mask_target%nlon)
    ierr = nf90_inq_dimid(ncid,'y',dimid)
    ierr = nf90_inquire_dimension(ncid,dimid,len=mask_target%nlat)
    WRITE(*,*) "Coupling dimensions ocean model grid: NX =" , mask_target%nlon, "  NY = ", mask_target%nlat
    ierr = nf90_close(ncid)

  END SUBROUTINE open_coupling_file
!
! *********************************************************************
  SUBROUTINE read_coupling_file(dn_couple, mask_src_mapped, ix_target, iy_target, mask_on_ocean)
! *********************************************************************
!
!   *** HD coupling file variables
!   *** HD grid 
!   ***   FMOU_HD_NEMO = Mask with HD mouth boxes that have an associated ocean model mouth point.
!   ***   INDEXX = x-Indices of nearest ocean model mouth points
!   ***   INDEXY = y-Indices of nearest ocean model mouth points
!   *** Ocean model grid 
!   *** FMOU_HD_ON_NEMO = Mask on which with HD mouth boxes are mapped

    CHARACTER (LEN=*), INTENT(in)   :: dn_couple              ! Coupling file name
    INTEGER, DIMENSION(:,:), INTENT(out) :: mask_src_mapped   ! Stores FMOU_HD_NEMO
    INTEGER, DIMENSION(:,:), INTENT(out) :: ix_target         ! Stores INDEXX
    INTEGER, DIMENSION(:,:), INTENT(out) :: iy_target         ! Stores INDEXY
!   Ocean model array
    INTEGER, DIMENSION(:,:), INTENT(out) :: mask_on_ocean     ! Stores FMOU_HD_ON_NEMO
!
    INTEGER :: ncid, varid, varid2, ierr  

    ierr = nf90_open(dn_couple, NF90_NOWRITE, ncid)

    ierr = nf90_inq_varid(ncid,'FMOU_HD_TO_NEMO',varid)
    ierr = nf90_get_var(ncid, varid, mask_src_mapped)
    ierr = nf90_inq_varid(ncid,'INDEXX',varid)
    ierr = nf90_get_var(ncid, varid, ix_target)
    ierr = nf90_inq_varid(ncid,'INDEXY',varid)
    ierr = nf90_get_var(ncid, varid, iy_target)

    WRITE(*,*) TRIM(dn_couple)," HD arrays read"

    ierr = nf90_inq_varid(ncid,'FMOU_HD_ON_NEMO',varid2)
    ierr = nf90_get_var(ncid, varid2, mask_on_ocean(:,:))
    ierr = nf90_close(ncid)
    WRITE(*,*) TRIM(dn_couple)," read and closed."

  END SUBROUTINE read_coupling_file
!
END MODULE mo_convert
