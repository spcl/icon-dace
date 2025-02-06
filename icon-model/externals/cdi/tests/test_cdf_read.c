#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "cdi.h"
#include "dmemory.h"

static void printAtts(int vlistID);

int
main(int argc, const char **argv)
{
  // todo: handle optional arguments here to increase test coverage
  const char *fname = (argc > 1) ? argv[1] : "test.nc";

  int countMissingValues = 1;

  int streamID = streamOpenRead(fname);
  if (streamID < 0)
    {
      fprintf(stderr, "Open failed for file %s: %s\n", fname, cdiStringError(streamID));
      return EXIT_FAILURE;
    }
  int vlistID = streamInqVlist(streamID);
  size_t nVars = (size_t) vlistNvars(vlistID);

  double *buf = NULL;
  size_t bufSize = 0;
  size_t allNmissSum = 0;

  printAtts(vlistID);

  for (int tsID = 0; streamInqTimestep(streamID, tsID); ++tsID)
    {
      for (size_t varID = 0; varID < nVars; ++varID)
        {
          size_t memSize = (size_t) vlistInqVarSize(vlistID, (int) varID) * sizeof(double);
          SizeType nmiss;
          if (memSize > bufSize)
            {
              double *temp = (double *) Realloc(buf, memSize);
              if (!temp)
                {
                  perror("read buffer reallocation failed");
                  return EXIT_FAILURE;
                }
              buf = temp;
            }
          streamReadVar(streamID, (int) varID, buf, &nmiss);
          allNmissSum += nmiss;
        }
    }
  if (countMissingValues) printf("missing values count = %zu\n", allNmissSum);
  free(buf);
  streamClose(streamID);
  return EXIT_SUCCESS;
}

static void
printAtts(int vlistID)
{
  size_t nVars = (size_t) vlistNvars(vlistID);
  static const char globDesc[] = "global", varDescPrefix[] = "variable ";
  size_t digitsPerInt = 9;
  size_t bufSize = digitsPerInt + sizeof(varDescPrefix);
  void *restrict buf = Malloc(bufSize);
  if (!buf)
    {
      perror("attribute buffer resize failed");
      exit(EXIT_FAILURE);
    }
  for (size_t varIdx = 0; varIdx <= nVars; ++varIdx)
    {
      int nAtts, attType, attLen;
      char attName[CDI_MAX_NAME + 1];
      int varID = (int) varIdx - 1;
      cdiInqNatts(vlistID, varID, &nAtts);
      for (size_t attIdx = 0; attIdx < (size_t) nAtts; attIdx++)
        {
          int rc = cdiInqAtt(vlistID, varID, (int) attIdx, attName, &attType, &attLen);
          {
            const char *varDesc
                = varIdx > 0
                      ? (snprintf((char *) buf, bufSize, "%s%d", varDescPrefix, vlistInqVarCode(vlistID, varID)), (char *) buf)
                      : globDesc;
            printf("%s attribute \"%s\", value: ", varDesc, attName);
          }
          if (attLen < 0)
            {
              puts("error retrieving value");
              continue;
            }
          size_t elemSize = 0;
          switch (attType)
            {
            case CDI_DATATYPE_TXT: elemSize = 1; break;
            case CDI_DATATYPE_FLT: elemSize = sizeof(double); break;
            case CDI_DATATYPE_INT: elemSize = sizeof(int); break;
            }

          size_t attSize = elemSize * ((size_t) attLen + 1);
          if (attSize > bufSize)
            {
              if (!(buf = Realloc(buf, attSize)))
                {
                  perror("attribute buffer resize failed");
                  exit(EXIT_FAILURE);
                }
            }

          switch (attType)
            {
            case CDI_DATATYPE_TXT:
              rc = cdiInqAttTxt(vlistID, (int) varID, attName, attLen, (char *) buf);
              if (rc == CDI_NOERR) printf("\"%.*s\"", attLen, (char *) buf);
              break;
            case CDI_DATATYPE_FLT:
              rc = cdiInqAttFlt(vlistID, (int) varID, attName, attLen + 1, (double *) buf);
              if (rc == CDI_NOERR && attLen)
                {
                  const double *restrict dp = (const double *) buf;
                  printf("%10g", dp[0]);
                  for (size_t i = 1; i < (size_t) attLen; ++i) printf(", %10g", dp[i]);
                }
              break;
            case CDI_DATATYPE_INT:
              rc = cdiInqAttInt(vlistID, (int) varID, attName, attLen + 1, (int *) buf);
              if (rc == CDI_NOERR && attLen)
                {
                  const int *restrict ip = (const int *) buf;
                  printf("%d", ip[0]);
                  for (size_t i = 1; i < (size_t) attLen; ++i) printf(", %d", ip[i]);
                }
              break;
            }

          if (rc == CDI_NOERR)
            {
              putchar('\n');
            }
          else
            {
              puts("error retrieving value");
            }
        }
    }
  free(buf);
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
