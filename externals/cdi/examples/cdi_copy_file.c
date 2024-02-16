#include <stdio.h>
#include <stdlib.h>

#include "cdi.h"

/* gcc -g -O2 -o cdi_copy_file cdi_copy_file.c -I../src ../src/libcdi.a -lm */

#define MAX_LEVEL 200

int
main(int argc, char *argv[])
{
  int linfo = 0;
  int vlistID1, vlistID2, varID, streamID1, streamID2;
  size_t nmiss;
  int nvars, status;
  double *vardata = NULL;
  const char *ifile, *ofile;

  if (argc == 1)
    {
      fprintf(stderr, "usage: %s ifile ofile\n", argv[0]);
      return (-1);
    }

  ifile = argv[1];
  ofile = argv[2];

  /* Open the input dataset */
  streamID1 = streamOpenRead(ifile);
  if (streamID1 < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID1));
      return (1);
    }

  /* Get the variable list of the dataset */
  vlistID1 = streamInqVlist(streamID1);
  nvars = vlistNvars(vlistID1);
  size_t gridsize = vlistGridsizeMax(vlistID1);

  vardata = (double *) malloc(gridsize * MAX_LEVEL * sizeof(double));

  /* Open the output dataset (GRIB format) */
  streamID2 = streamOpenWrite(ofile, CDI_FILETYPE_GRB);
  if (streamID2 < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID2));
      return (1);
    }

  vlistID2 = vlistDuplicate(vlistID1);

  streamDefVlist(streamID2, vlistID2);

  /* Loop over all time steps */
  int tsID = 0;
  while ((status = streamInqTimestep(streamID1, tsID)))
    {
      /* Define the output time step */
      streamDefTimestep(streamID2, tsID);

      for (varID = 0; varID < nvars; ++varID)
        {
          /* Read variable  */
          streamReadVar(streamID1, varID, vardata, &nmiss);

          if (linfo)
            {
              int i, k, nlevel;
              nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
              gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
              for (k = 0; k < nlevel; ++k)
                {
                  double fmin = 1.e30, fmax = -1.e30, fmean = 0;
                  for (i = 0; i < gridsize; ++i)
                    {
                      if (vardata[k * gridsize + i] < fmin) fmin = vardata[k * gridsize + i];
                      if (vardata[k * gridsize + i] > fmax) fmax = vardata[k * gridsize + i];
                      fmean += vardata[k * gridsize + i];
                    }
                  printf("%3d %3d %3d %9g %9g %9g\n", tsID, varID, k, fmin, fmean / gridsize, fmax);
                }
            }

          /* Write variable */
          streamWriteVar(streamID2, varID, vardata, nmiss);
        }

      tsID++;
    }

  /* Close the streams */
  streamClose(streamID1);
  streamClose(streamID2);

  free(vardata);

  return 0;
}
