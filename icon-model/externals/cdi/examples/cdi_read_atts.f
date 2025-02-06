      PROGRAM CDIREAD

      IMPLICIT NONE

      INCLUDE 'cdi.inc'

      INTEGER nlon, nlat, nlev, nts
      PARAMETER (nlon = 12)   ! Number of longitudes
      PARAMETER (nlat =  6)   ! Number of latitudes
      PARAMETER (nlev =  5)   ! Number of levels
      PARAMETER (nts  =  3)   ! Number of time steps

      CHARACTER(80) varname
      CHARACTER(80) attname
      CHARACTER(80) atttxt
      INTEGER ia, natts, atttype, attlen
      INTEGER inst
      INTEGER gridID, zaxisID1, zaxisID2, taxisID
      INTEGER vlistID, varID1, varID2, streamID, tsID
      INTEGER status, vdate, vtime
      INTEGER nmiss
      REAL*8 var1(nlon*nlat), var2(nlon*nlat*nlev)

!     Open the dataset
      streamID = streamOpenRead("example2.nc")
      IF ( streamID < 0 ) THEN
         WRITE(0,*) cdiStringError(streamID)
         STOP
      END IF

!     Get the variable list of the dataset
      vlistID = streamInqVlist(streamID)
      varname(1:80) = "       "
      CALL vlistInqVarName(vlistID, 0, varname)
      WRITE(*,*) 'varname : ', varname

!     Set the variable IDs
      varID1 = 0
      varID2 = 1

!     Get the Time axis from the variable list
      taxisID = vlistInqTaxis(vlistID)

      status = vlistInqNatts(vlistID, -1, natts);
      WRITE(0,*) 'natts: ', natts
      attname(1:80) = " "

      DO ia = 1, natts

        status = vlistInqAtt(vlistID, -1, ia-1, attname, atttype,
     &  attlen)
        IF ( atttype == CDI_DATATYPE_TXT ) THEN
	  status = vlistInqAttTxt(vlistID, -1, attname, 80,
     &         atttxt)
          WRITE(0,*) attname(1:10), attlen, atttxt(1:attlen)
!	  atttxt[attlen] = 0;
!	  fprintf(fp, "  %s=\"%s\"\n", attname, atttxt);
        ENDIF
      END DO

!     Loop over the number of time steps
      DO tsID = 0, nts-1
!        Inquire the time step
         status = streamInqTimestep(streamID, tsID)

!        Get the verification date and time
         vdate = taxisInqVdate(taxisID)
         vtime = taxisInqVtime(taxisID)

!        Read var1 and var2
         CALL streamReadVar(streamID, varID1, var1, nmiss)
         CALL streamReadVar(streamID, varID2, var2, nmiss)
      END DO

!     Close the input stream
      CALL streamClose(streamID)

!      CALL cdiReset()

      END
