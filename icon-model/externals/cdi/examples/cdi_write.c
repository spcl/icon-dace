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
  SizeType numMissVals = 0;
  double lons[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
  double lats[] = { -75, -45, -15, 15, 45, 75 };
  double levels[] = { 101300, 92500, 85000, 50000, 20000 };
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Create a regular lon/lat grid
  int gridID = gridCreate(GRID_LONLAT, nlon * nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);

  // Create a surface level Z-axis
  int zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a pressure level Z-axis
  int zaxisID2 = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID2, levels);

  // Create a variable list
  int vlistID = vlistCreate();

  // Define the variables
  int varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TIME_VARYING);
  int varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARYING);

  // Define the variable names
  vlistDefVarName(vlistID, varID1, "varname1");
  vlistDefVarName(vlistID, varID2, "varname2");

  // Create a Time axis
  int taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset in netCDF format
  int streamID = streamOpenWrite("example.nc", CDI_FILETYPE_NC);
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return 1;
    }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  // Loop over the number of time steps
  for (int tsID = 0; tsID < nts; tsID++)
    {
      // Set the verification date to 1985-01-01 + tsID
      taxisDefVdate(taxisID, 19850101 + tsID);
      // Set the verification time to 12:00:00
      taxisDefVtime(taxisID, 120000);
      // Define the time step
      streamDefTimestep(streamID, tsID);

      // Init var1 and var2
      for (size_t i = 0; i < nlon * nlat; i++) var1[i] = 1.1;
      for (size_t i = 0; i < nlon * nlat * nlev; i++) var2[i] = 2.2;

      // Write var1 and var2
      streamWriteVar(streamID, varID1, var1, numMissVals);
      streamWriteVar(streamID, varID2, var2, numMissVals);
    }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID1);
  zaxisDestroy(zaxisID2);
  gridDestroy(gridID);

  return 0;
}
