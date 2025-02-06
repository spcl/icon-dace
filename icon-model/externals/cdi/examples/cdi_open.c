#include <stdio.h>
#include <stdlib.h>
#include "cdi.h"

#define NOPEN 2048
#define FILENAME "example.grb"

int
main(void)
{
  int iopen;
  int idel;
  int taxisID, vlistID, streamID;
  int streamIDs[NOPEN];

  for (iopen = 0; iopen < NOPEN; ++iopen)
    {
      fprintf(stderr, "sequential open %d\n", iopen + 1);
      /* Open the dataset */
      streamID = streamOpenRead(FILENAME);
      if (streamID < 0)
        {
          fprintf(stderr, "%s\n", cdiStringError(streamID));
          return (1);
        }

      /* Get the variable list of the dataset */
      vlistID = streamInqVlist(streamID);

      /* Get the Time axis from the variable list */
      taxisID = vlistInqTaxis(vlistID);

      /* Close the input stream */
      streamClose(streamID);
    }

  for (iopen = 0; iopen < NOPEN; ++iopen) streamIDs[iopen] = -1;

  for (iopen = 0; iopen < NOPEN; ++iopen)
    {
      fprintf(stderr, "simultaneous open %d\n", iopen + 1);
      /* Open the dataset */
      streamID = streamOpenRead(FILENAME);
      if (streamID < 0)
        {
          fprintf(stderr, "%s\n", cdiStringError(streamID));
          return (1);
        }

      /* Get the variable list of the dataset */
      vlistID = streamInqVlist(streamID);

      /* Get the Time axis from the variable list */
      taxisID = vlistInqTaxis(vlistID);

      streamIDs[iopen] = streamID;

      if (iopen % 2)
        {
          idel = (int) ((iopen + 1) * (rand() / (RAND_MAX + 1.0)));
          if (streamIDs[idel] >= 0)
            {
              streamClose(streamIDs[idel]);
              streamIDs[idel] = -1;
              fprintf(stderr, "randomly closed %d\n", idel);
            }
        }
    }

  for (iopen = 0; iopen < NOPEN; ++iopen)
    {
      /* Close the input stream */
      if (streamIDs[iopen] >= 0) streamClose(streamIDs[iopen]);
    }

  for (iopen = 0; iopen < NOPEN; ++iopen)
    {
      fprintf(stderr, "simultaneous open %d\n", iopen + 1);
      /* Open the dataset */
      streamID = streamOpenRead(FILENAME);
      if (streamID < 0)
        {
          fprintf(stderr, "%s\n", cdiStringError(streamID));
          return (1);
        }

      /* Get the variable list of the dataset */
      vlistID = streamInqVlist(streamID);

      /* Get the Time axis from the variable list */
      taxisID = vlistInqTaxis(vlistID);

      streamIDs[iopen] = streamID;
    }

  for (iopen = 0; iopen < NOPEN; ++iopen)
    {
      /* Close the input stream */
      if (streamIDs[iopen] >= 0) streamClose(streamIDs[iopen]);
    }

  return 0;
}
