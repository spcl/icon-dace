#include <stdio.h>
#include "cdi.h"

int
main(void)
{
  enum
  {
    nlon = 12,  // Number of longitudes
    nlat = 6,   // Number of latitudes
    nlev = 5,   // Number of levels
    nts = 3,    // Number of time steps
  };
  SizeType numMissVals;
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Open the dataset
  int streamID = streamOpenRead("example.nc");
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return 1;
    }

  // Get the variable list of the dataset
  int vlistID = streamInqVlist(streamID);

  // Set the variable IDs
  int varID1 = 0;
  int varID2 = 1;

  // Get the Time axis from the variable list
  int taxisID = vlistInqTaxis(vlistID);

  // Loop over the number of time steps
  for (int tsID = 0; tsID < nts; tsID++)
    {
      // Inquire the time step
      streamInqTimestep(streamID, tsID);

      // Get the verification date and time
      int vdate = taxisInqVdate(taxisID);
      int vtime = taxisInqVtime(taxisID);
      printf("read timestep %d:  date=%d  time=%d\n", tsID + 1, vdate, vtime);

      // Read var1 and var2
      streamReadVar(streamID, varID1, var1, &numMissVals);
      streamReadVar(streamID, varID2, var2, &numMissVals);
    }

  // Close the input stream
  streamClose(streamID);

  return 0;
}
