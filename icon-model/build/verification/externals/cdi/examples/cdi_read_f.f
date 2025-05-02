      PROGRAM CDIREAD

      IMPLICIT NONE

      INCLUDE 'cdi.inc'

      INTEGER nlon, nlat
      INTEGER nlev, nts
      PARAMETER (nlon = 12)   ! Number of longitudes
      PARAMETER (nlat =  6)   ! Number of latitudes
      PARAMETER (nlev =  5)   ! Number of levels
      PARAMETER (nts  =  3)   ! Number of time steps

      INTEGER gridID, zaxisID1, zaxisID2, taxisID
      INTEGER vlistID, varID1, varID2, streamID, tsID
      INTEGER status, vdate, vtime
      INTEGER nmiss
      REAL*8 var1(nlon*nlat), var2(nlon*nlat*nlev)

!     Open the dataset
      streamID = streamOpenRead("example.nc")
      IF ( streamID < 0 ) THEN
         WRITE(0,*) cdiStringError(streamID)
         STOP
      END IF

!     Get the variable list of the dataset
      vlistID = streamInqVlist(streamID)

!     Set the variable IDs
      varID1 = 0
      varID2 = 1

!     Get the Time axis from the variable list
      taxisID = vlistInqTaxis(vlistID)

!     Loop over the number of time steps
      DO tsID = 0, nts-1
!        Inquire the time step
         status = streamInqTimestep(streamID, tsID)

!        Get the verification date and time
         vdate = taxisInqVdate(taxisID)
         vtime = taxisInqVtime(taxisID)
         WRITE(0, *)"read timestep ",tsID+1,"date=",vdate,"time=",vtime


!        Read var1 and var2
         CALL streamReadVar(streamID, varID1, var1, nmiss)
         CALL streamReadVar(streamID, varID2, var2, nmiss)
      END DO

!     Close the input stream
      CALL streamClose(streamID)

      END
