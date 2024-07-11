#include <stdio.h>
#include <stdlib.h>
#include "cdi.h"

int
main(void)
{
  // Open the input dataset
  int streamID1 = streamOpenRead("example.nc");
  if (streamID1 < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID1));
      return 1;
    }

  // Get the variable list of the dataset
  int vlistID1 = streamInqVlist(streamID1);

  int numVars = vlistNvars(vlistID1);

  int varDataSize = 0;
  double *varData = NULL;
  for (int varID = 0; varID < numVars; ++varID)
    {
      int varSize = vlistInqVarSize(vlistID1, varID);
      varDataSize = varSize > varDataSize ? varSize : varDataSize;
    }
  varData = malloc((size_t) varDataSize * sizeof(double));
  if (!varData)
    {
      perror("cannot allocate temporary copying buffer");
      return EXIT_FAILURE;
    }

  // Open the output dataset (GRIB format)
  int streamID2 = streamOpenWrite("example.grb", CDI_FILETYPE_GRB);
  if (streamID2 < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID2));
      return EXIT_FAILURE;
    }

  int vlistID2 = vlistDuplicate(vlistID1);

  streamDefVlist(streamID2, vlistID2);

  // Loop over the input time steps
  int tsID = 0;
  while (streamInqTimestep(streamID1, tsID))
    {
      // Define the output time step
      streamDefTimestep(streamID2, tsID);

      for (int varID = 0; varID < numVars; ++varID)
        {
          SizeType numMissVals;
          // Read var
          streamReadVar(streamID1, varID, varData, &numMissVals);
          // Write var
          streamWriteVar(streamID2, varID, varData, numMissVals);
        }
      ++tsID;
    }

  // Close the streams
  streamClose(streamID1);
  streamClose(streamID2);

  return EXIT_SUCCESS;
}
