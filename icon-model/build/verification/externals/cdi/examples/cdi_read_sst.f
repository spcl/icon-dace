      PROGRAM CDIREADSST

      IMPLICIT NONE

      INCLUDE 'cdi.inc'

      INTEGER nvals
      PARAMETER (nvals = 20480)

      INTEGER gridID, taxisID
      INTEGER vlistID, varID1, streamID, tsID
      INTEGER status, vdate, vtime
      INTEGER nmiss
      REAL*8 sst(nvals)

!      CALL cdiDebug(1)
!     Open the dataset
      streamID = streamOpenRead
     & ("/Users/m214003/data/icon_amip2sst_1870-2010.nc")
      IF ( streamID < 0 ) THEN
         WRITE(0,*) cdiStringError(streamID)
         STOP
      END IF

!     Get the variable list of the dataset
      vlistID = streamInqVlist(streamID)

!     Set the variable IDs
      varID1 = 0

!     Get the Time axis from the variable list
      taxisID = vlistInqTaxis(vlistID)

!     Loop over the first 10 time steps
      DO tsID = 0, 10
!        Inquire the time step
         status = streamInqTimestep(streamID, tsID)

!        Get the verification date and time
         vdate = taxisInqVdate(taxisID)
         vtime = taxisInqVtime(taxisID)

!        Read sst
         CALL streamReadVarSlice(streamID, varID1, 0, sst, nmiss)

         WRITE(*,*) vdate, minval(sst), maxval(sst)
      END DO

!     Close the input stream
      CALL streamClose(streamID)

      END
