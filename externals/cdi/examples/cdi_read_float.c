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
  float var1[nlon * nlat];
  float var2[nlon * nlat * nlev];

  // Open the dataset
  streamID = streamOpenRead("example.grb");
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return (1);
    }

  // Get the variable list of the dataset
  vlistID = streamInqVlist(streamID);

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
      streamReadVarF(streamID, varID1, var1, &nmiss);
      streamReadVarF(streamID, varID2, var2, &nmiss);
      printf("tsID %d %d %d %g %g\n", tsID, vdate, vtime, var1[0], var2[0]);
    }

  // Close the input stream
  streamClose(streamID);

  return 0;
}
