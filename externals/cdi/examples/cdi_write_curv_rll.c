#include <stdio.h>
#include "cdi.h"

int
main(void)
{
  const int nlon = 12;  // Number of longitudes
  const int nlat = 6;   // Number of latitudes
  const int nlev = 5;   // Number of levels
  const int nts = 3;    // Number of time steps
  size_t nmiss = 0;
  double levs[] = { 101300, 92500, 85000, 50000, 20000 };
  double var1[nlon * nlat];
  double var2[nlon * nlat * nlev];

  // Create a curvilinear grid
  int gridID = gridCreate(GRID_CURVILINEAR, nlon * nlat);
  {
    double lons[] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
    double lats[] = { -75, -45, -15, 15, 45, 75 };

    double lons2D[nlon * nlat], lats2D[nlon * nlat];
    for (int j = 0; j < nlat; ++j)
      for (int i = 0; i < nlon; ++i)
        {
          lons2D[j * nlon + i] = lons[i];
          lats2D[j * nlon + i] = lats[j];
        }

    gridDefXsize(gridID, nlon);
    gridDefYsize(gridID, nlat);
    gridDefXvals(gridID, lons2D);
    gridDefYvals(gridID, lats2D);

    cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, "lon");
    cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, "lat");
  }

  // Create a rotated lon/lat grid
  {
    int projID = gridCreate(GRID_PROJECTION, nlon * nlat);
    double rlons[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 };
    double rlats[] = { 5, 10, 15, 20, 25, 30 };

    gridDefXsize(projID, nlon);
    gridDefYsize(projID, nlat);
    gridDefXvals(projID, rlons);
    gridDefYvals(projID, rlats);

    cdiDefKeyString(projID, CDI_XAXIS, CDI_KEY_NAME, "rlon");
    cdiDefKeyString(projID, CDI_YAXIS, CDI_KEY_NAME, "rlat");

    gridDefParamRLL(projID, -162., 39.25, 0);

    gridDefProj(gridID, projID);
  }

  // Create a surface level Z-axis
  int zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a pressure level Z-axis
  int zaxisID2 = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID2, levs);

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
  int streamID = streamOpenWrite("example.grb", CDI_FILETYPE_GRB2);
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
      streamWriteVar(streamID, varID1, var1, nmiss);
      streamWriteVar(streamID, varID2, var2, nmiss);
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
