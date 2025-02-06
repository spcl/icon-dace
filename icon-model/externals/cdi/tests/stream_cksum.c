#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>

#include "cdi.h"
#include "cksum.h"
#include "stream_cksum.h"
#include "dmemory.h"

const char *zaxisNamePtr(int zaxistype);

struct cksum_table *
cksum_stream(const char *fname, size_t *table_len)
{
  int nvars;
  uint32_t *checksum_state = NULL;
  enum directionZ
  {
    DIRECTION_DOWN = 1,
    DIRECTION_NONE,
    DIRECTION_UP
  };
  struct var_desc_t
  {
    int x, y, z;
    enum directionZ zDirection;
    int code;
    size_t chars;
  } *varDesc = NULL;
  size_t var_size_max_chars = 0;
  double *buf = NULL;
  struct cksum_table *file_vars = NULL;

  do
    {
      // Open the dataset
      int streamID = streamOpenRead(fname);
      if (streamID < 0)
        {
          fprintf(stderr, "Cannot open data input file %s: %s\n", fname, cdiStringError(streamID));
          nvars = -1;
          break;
        }

      // Get the variable list of the dataset
      int vlistID = streamInqVlist(streamID);
      int fileType = streamInqFiletype(streamID);
      bool isLegacyFile = fileType == CDI_FILETYPE_SRV || fileType == CDI_FILETYPE_EXT;
      nvars = vlistNvars(vlistID);
      int ngrids = vlistNgrids(vlistID);
      int nzaxis = vlistNzaxis(vlistID);
      if (nzaxis < 0 || ngrids < 0)
        {
          fprintf(stderr, "Error in grid/zaxis count query %d:%d\n", ngrids, nzaxis);
          nvars = -1;
          break;
        }
      checksum_state = (uint32_t *) Calloc((size_t) nvars, sizeof(checksum_state[0]));
      varDesc = (struct var_desc_t *) Malloc((size_t) nvars * sizeof(varDesc[0]));

      for (int varIdx = 0; varIdx < nvars; ++varIdx)
        {
          int grid = vlistInqVarGrid(vlistID, varIdx), gridType, varCode = vlistInqVarCode(vlistID, varIdx);
          varDesc[varIdx].code = varCode;
          int zaxisID = vlistInqVarZaxis(vlistID, varIdx);
          if (grid == CDI_UNDEFID || zaxisID == CDI_UNDEFID)
            {
              fputs("error in axis/grid inquiry\n", stderr);
              nvars = -1;
              break;
            }
          int zSize;
          if ((zSize = varDesc[varIdx].z = zaxisInqSize(zaxisID)) <= 0)
            {
              fputs("invalid Z-axis found\n", stderr);
              nvars = -1;
              break;
            }
          if (isLegacyFile)
            varDesc[varIdx].zDirection = DIRECTION_NONE;
          else if (zSize > 1)
            {
              double lev[2];
              for (int levIdx = 0; levIdx < 2; ++levIdx) lev[levIdx] = zaxisInqLevel(zaxisID, levIdx);
              int zaxistype = zaxisInqType(zaxisID);
              switch (zaxistype)
                {
                case ZAXIS_PRESSURE:
                  if (lev[0] < lev[1])
                    varDesc[varIdx].zDirection = DIRECTION_DOWN;
                  else if (lev[1] < lev[0])
                    varDesc[varIdx].zDirection = DIRECTION_UP;
                  else
                    {
                      fprintf(stderr,
                              "unexpected level ordering on z-Axis for variable"
                              " code %d found: lev[0]=%g, lev[1]=%g\n",
                              varCode, lev[0], lev[1]);
                      nvars = -1;
                    }
                  break;
                default:
                  fprintf(stderr,
                          "unexpected type of z-Axis for variable"
                          " code %d found: %s\n",
                          varCode, zaxisNamePtr(zaxistype));
                  nvars = -1;
                }
            }
          else
            varDesc[varIdx].zDirection = DIRECTION_NONE;
          if (nvars == -1) break;
          if ((gridType = gridInqType(grid)) != GRID_LONLAT && gridType != GRID_GENERIC && gridType != GRID_CURVILINEAR)
            {
              fprintf(stderr, "unexpected non-lonlat grid found: %d\n", gridType);
              nvars = -1;
              break;
            }
          if ((varDesc[varIdx].x = gridInqXsize(grid)) < 0)
            {
              fprintf(stderr, "invalid X-size found: %d\n", varDesc[varIdx].x);
              nvars = -1;
              break;
            }
          if (varDesc[varIdx].x == 0) varDesc[varIdx].x = 1;
          if ((varDesc[varIdx].y = gridInqYsize(grid)) < 0)
            {
              fprintf(stderr, "invalid Y-size found: %d\n", varDesc[varIdx].y);
              nvars = -1;
              break;
            }
          if (varDesc[varIdx].y == 0) varDesc[varIdx].y = 1;
          varDesc[varIdx].chars
              = (size_t) varDesc[varIdx].x * (size_t) varDesc[varIdx].y * (size_t) varDesc[varIdx].z * sizeof(buf[0]);
          if (var_size_max_chars < varDesc[varIdx].chars) var_size_max_chars = varDesc[varIdx].chars;
        }
      buf = (double *) Malloc(var_size_max_chars);

      if (nvars == -1) break;

      // Get the Time axis from the variable list
      int taxisID = vlistInqTaxis(vlistID);

      int tsID = 0;
      // Inquire the time step
      while (streamInqTimestep(streamID, tsID))
        {
          // Get the verification date and time
          int vdatetime[2] = { taxisInqVtime(taxisID), (int) taxisInqVdate(taxisID) };
          // Read var1 and var2
          for (int varIdx = 0; varIdx < nvars; ++varIdx)
            {
              SizeType nmiss;
              streamReadVar(streamID, varIdx, buf, &nmiss);
              memcrc_r(checksum_state + varIdx, (const unsigned char *) vdatetime, sizeof(vdatetime));
              if (varDesc[varIdx].zDirection == DIRECTION_UP || varDesc[varIdx].zDirection == DIRECTION_NONE)
                memcrc_r(checksum_state + varIdx, (const unsigned char *) buf, varDesc[varIdx].chars);
              else
                {
                  size_t nlev = (size_t) varDesc[varIdx].z,
                         charsPerLev = (size_t) varDesc[varIdx].x * (size_t) varDesc[varIdx].y * sizeof(buf[0]);
                  for (size_t lev = 0; lev < nlev; ++lev)
                    memcrc_r(checksum_state + varIdx, (const unsigned char *) buf + (nlev - lev - 1) * charsPerLev, charsPerLev);
                }
            }
          ++tsID;
        }

      file_vars = (struct cksum_table *) Malloc((size_t) nvars * sizeof(file_vars[0]));
      for (int varIdx = 0; varIdx < nvars; ++varIdx)
        {
          file_vars[varIdx].code = varDesc[varIdx].code;
          file_vars[varIdx].cksum
              = memcrc_finish(checksum_state + varIdx, (off_t) ((varDesc[varIdx].chars + sizeof(int) * 2) * (size_t) tsID));
        }
      // Close the input stream
      streamClose(streamID);
    }
  while (0);

  // free resources
  free(checksum_state);
  free(varDesc);
  free(buf);
  *table_len = (size_t) nvars;

  return file_vars;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
