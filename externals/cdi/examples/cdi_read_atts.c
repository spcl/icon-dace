#include <stdio.h>
#include "cdi.h"

int nlon = 12;  // Number of longitudes
int nlat = 6;   // Number of latitudes
int nlev = 5;   // Number of levels
int nts = 3;    // Number of time steps

int
main(void)
{
  int taxisID, vlistID, varID1, varID2, streamID;
  size_t nmiss;
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Open the dataset
  streamID = streamOpenRead("example_att.nc");
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return (1);
    }

  // Get the variable list of the dataset
  vlistID = streamInqVlist(streamID);

  int status;
  int natts;
  status = vlistInqNatts(vlistID, -1, &natts);
  printf("vlistInqNatts status %d\n", status);
  int type, len;
  char attname[256];
  for (int i = 0; i < natts; ++i)
    {
      status = vlistInqAtt(vlistID, -1, i, attname, &type, &len);
      printf("vlistInqAtt status = %d  name = %s  len = %d\n", status, attname, len);
    }

  int ival;
  status = vlistInqAttInt(vlistID, -1, "iii", 1, &ival);
  printf("vlistInqAttInt status = %d  name = %s  val = %d\n", status, "iii", ival);
  status = vlistInqAttInt(vlistID, -1, "kkk", 1, &ival);
  printf("vlistInqAttInt status = %d  name = %s  val = %d\n", status, "kkk", ival);
  status = vlistInqAttInt(vlistID, -1, "jjj", 1, &ival);
  printf("vlistInqAttInt status = %d  name = %s  val = %d\n", status, "jjj", ival);

  // Set the variable IDs
  varID1 = 0;
  varID2 = 1;

  // Get the Time axis from the variable list
  taxisID = vlistInqTaxis(vlistID);

  // Loop over the number of time steps
  for (int tsID = 0; tsID < nts; tsID++)
    {
      // Inquire the time step
      streamInqTimestep(streamID, tsID);

      // Get the verification date and time
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);

      // Read var1 and var2
      streamReadVar(streamID, varID1, var1, &nmiss);
      streamReadVar(streamID, varID2, var2, &nmiss);
    }

  // Close the input stream
  streamClose(streamID);

  return 0;
}
