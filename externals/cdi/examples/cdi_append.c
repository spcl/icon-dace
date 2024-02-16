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
  SizeType nmiss = 0;
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Append a dataset (created by cdi_write)
  int streamID = streamOpenAppend("example.nc");
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return 1;
    }

  // Get the variable list of the dataset
  int vlistID = streamInqVlist(streamID);

  int numVars = vlistNvars(vlistID);
  if (numVars != 2)
    {
      fprintf(stderr, "Unexpected number of variables: %d\n", numVars);
      return 1;
    }

  int varID1 = 0;
  int varID2 = 1;

  SizeType varSize1 = gridInqSize(vlistInqVarGrid(vlistID, varID1)) * zaxisInqSize(vlistInqVarZaxis(vlistID, varID1));
  if (varSize1 != nlon * nlat)
    {
      fprintf(stderr, "Unexpected size of variable 1: %ld\n", (long) varSize1);
      return 1;
    }

  SizeType varSize2 = gridInqSize(vlistInqVarGrid(vlistID, varID2)) * zaxisInqSize(vlistInqVarZaxis(vlistID, varID2));
  if (varSize2 != nlon * nlat * nlev)
    {
      fprintf(stderr, "Unexpected size of variable 2: %ld\n", (long) varSize2);
      return 1;
    }

  int taxisID = vlistInqTaxis(vlistID);

  int nsteps = vlistNtsteps(vlistID);
  printf("Number of timesteps: %d\n", nsteps);

  // Loop over the number of time steps
  for (int tsID = 0; tsID < nts; tsID++)
    {
      // Set the verification date to 1985-01-01 + tsID
      taxisDefVdate(taxisID, 19850101 + nsteps + tsID);
      // Set the verification time to 12:00:00
      taxisDefVtime(taxisID, 120000);
      // Define the time step
      streamDefTimestep(streamID, nsteps + tsID);

      // Init var1 and var2
      for (size_t i = 0; i < nlon * nlat; i++) var1[i] = 1.1;
      for (size_t i = 0; i < nlon * nlat * nlev; i++) var2[i] = 2.2;

      // Write var1 and var2
      streamWriteVar(streamID, varID1, var1, nmiss);
      streamWriteVar(streamID, varID2, var2, nmiss);
    }

  // Close the output stream
  streamClose(streamID);

  return 0;
}
