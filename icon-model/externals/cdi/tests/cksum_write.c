#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cdi.h"
#include "cdi_uuid.h"
#include "dmemory.h"

#include "cksum.h"
#include "simple_model_helper.h"

static int
parse_intarg(const char msg[])
{
  char *end;
  long temp = strtol(optarg, &end, 0);
  if ((errno == ERANGE && (temp == LONG_MAX || temp == LONG_MIN)) || (errno != 0 && temp == 0))
    {
      perror(msg);
      exit(EXIT_FAILURE);
    }
  if (temp > INT_MAX || temp < INT_MIN)
    {
      fprintf(stderr, "range error: %ld\n", temp);
      exit(EXIT_FAILURE);
    }
  return (int) temp;
}

/* If we're not using GNU C, elide __attribute__ */
#if !defined __GNUC__ && !defined __attribute__
#define __attribute__(x) /*NOTHING*/
#endif

static inline void
check_positive(int v, const char *msg)
{
  if (v < 1)
    {
      fprintf(stderr, "error: number of %s must be positive!\n", msg);
      exit(EXIT_FAILURE);
    }
}

#ifdef TEST_CHUNK_WRITE
static void
get_chunk(double *chunkBuf, double *var, int varShape[3], int chunk[3][2])
{
  size_t ofs = 0;
  size_t start_k = (size_t) chunk[2][0], start_j = (size_t) chunk[1][0], start_i = (size_t) chunk[0][0];
  size_t size_k = (size_t) chunk[2][1] - (size_t) chunk[2][0] + 1, size_j = (size_t) chunk[1][1] - (size_t) chunk[1][0] + 1,
         size_i = (size_t) chunk[0][1] - (size_t) chunk[0][0] + 1;
  size_t stride_k = (size_t) varShape[0] * (size_t) varShape[1], stride_j = (size_t) varShape[0];
  for (size_t k = 0; k < size_k; ++k)
    for (size_t j = 0; j < size_j; ++j)
      for (size_t i = 0; i < size_i; ++i)
        chunkBuf[ofs++] = var[(k + start_k) * stride_k + (j + start_j) * stride_j + (i + start_i)];
}
#endif

#ifndef TEST_CHUNK_WRITE
static const struct
{
  char suffix[5];
  int type, defaultDT, defaultGrid;
} suffix2type[] = {
  { "nc", CDI_FILETYPE_NC, CDI_DATATYPE_FLT64, GRID_LONLAT },
  { "grb", CDI_FILETYPE_GRB, CDI_DATATYPE_PACK24, GRID_LONLAT },
  { "grb2", CDI_FILETYPE_GRB2, CDI_DATATYPE_PACK24, GRID_LONLAT },
  { "nc2", CDI_FILETYPE_NC2, CDI_DATATYPE_FLT64, GRID_LONLAT },
  { "nc4", CDI_FILETYPE_NC4, CDI_DATATYPE_FLT64, GRID_LONLAT },
  {
      "ext",
      CDI_FILETYPE_EXT,
      CDI_DATATYPE_FLT64,
      GRID_GENERIC,
  },
  {
      "srv",
      CDI_FILETYPE_SRV,
      CDI_DATATYPE_FLT64,
      GRID_GENERIC,
  },
  { "ieg", CDI_FILETYPE_IEG, CDI_DATATYPE_FLT64, GRID_LONLAT },
};
#endif
enum
{
  nvars = 2,
};

static const int varCodes[nvars] = { 42, 55 };

int
main(int argc, char *argv[])
{
  int gridID, zaxisID[nvars], taxisID;
  int vlistID, varID[nvars], streamID;
  int nlon = 12,  //!< Number of longitudes
      nlat = 6,   //!< Number of latitudes
      nlev = 5,   //!< Number of levels
      nts = 3;    //!< Number of time steps
  enum
  {
    nmiss = 0
  };
  double *restrict var[nvars], mscale, mrscale;
  size_t varSize[nvars];
  const char *varName[nvars] = { "varname1", "varname2" };
#ifndef TEST_CHUNK_WRITE
  const char *suffix = "grb", *prefix = "example";
  int filetype = CDI_FILETYPE_GRB, datatype = CDI_DATATYPE_PACK24;
#else
  const char *suffix = "nc", *prefix = "example";
  int filetype = CDI_FILETYPE_NC, datatype = CDI_DATATYPE_FLT64;
#endif
  {
    int opt;
    while ((opt = getopt(argc, argv,
#ifndef TEST_CHUNK_WRITE
                         "f:"
#endif
                         "b:m:n:o:t:"))
           != -1)
      switch (opt)
        {
#ifndef TEST_CHUNK_WRITE
        case 'f':
          {
            int found = 0;
            for (size_t i = 0; i < sizeof(suffix2type) / sizeof(suffix2type[0]); ++i)
              if (!strcmp(optarg, suffix2type[i].suffix))
                {
                  found = 1;
                  filetype = suffix2type[i].type;
                  suffix = suffix2type[i].suffix;
                  datatype = suffix2type[i].defaultDT;
                  break;
                }
            if (!found)
              {
                fprintf(stderr, "Unsupported format requested: %s\n", optarg);
                exit(EXIT_FAILURE);
              }
          }
          break;
#endif
        case 'b': prefix = optarg; break;
        case 'm': nlon = parse_intarg("error parsing number of longitudes"); check_positive(nlon, "longitudes");
#ifdef TEST_CHUNK_WRITE
          if (nlon < 2)
            {
              fputs("number of longitudes must be larger 1 for "
                    "chunk write test\n",
                    stderr);
              exit(EXIT_FAILURE);
            }
#endif
          break;
        case 'n': check_positive(nlat = parse_intarg("error parsing number of latitudes"), "latitudes"); break;
        case 'o': check_positive(nlev = parse_intarg("error parsing number of levels"), "levels"); break;
        case 't': check_positive(nts = parse_intarg("error parsing number of timesteps"), "timesteps"); break;
        default: /* '?' */ fprintf(stderr, "Usage: %s [-m nlon] [-n nlat] [-o nlev] [-t nts]\n", argv[0]); exit(EXIT_FAILURE);
        }
  }

  varSize[0] = (size_t) nlon * (size_t) nlat;
  varSize[1] = (size_t) nlon * (size_t) nlat * (size_t) nlev;

  // Create a regular lon/lat grid
  gridID = createGlobalLatLonGrid(nlon, nlat);

  // Create a surface level Z-axis
  zaxisID[0] = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a pressure level Z-axis
  zaxisID[1] = zaxisCreate(ZAXIS_PRESSURE, nlev);
  {
    double *levs = (double *) Malloc((size_t) nlev * sizeof(levs[0]));
    for (size_t i = 0; i < (size_t) nlev; ++i) levs[i] = 101300 - floor(3940.3 * expm1(2.3579 * (double) i / (nlev - 1)));
    zaxisDefLevels(zaxisID[1], levs);
    free(levs);
  }

  /* add uuids to zaxis and grid */
  {
    unsigned char uuid[CDI_UUID_SIZE];

    static char gridUUIDTxt[] = "107d7a5b-348c-4d1a-90a9-d745914f2fb6";

    cdiStr2UUID(gridUUIDTxt, uuid);
    gridDefUUID(gridID, uuid);

    static char zaxisUUIDTxt[2][37] = { { "d157f399-5496-4097-a3d8-437a6dda6311" }, { "6f784a65-bce8-48c9-afa4-4c40130709c7" } };

    for (int i = 0; i < 2; ++i)
      {
        cdiStr2UUID(zaxisUUIDTxt[i], uuid);
        cdiDefKeyBytes(zaxisID[i], CDI_GLOBAL, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);
      }
  }

  // Create a Time axis
  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Create a variable list
  vlistID = vlistCreate();

  for (size_t i = 0; i < nvars; ++i)
    {
      // Define the variables
      varID[i] = vlistDefVar(vlistID, gridID, zaxisID[i], TIME_VARIABLE);
      // Define the variable names,
      vlistDefVarName(vlistID, varID[i], varName[i]);
      // the codes
      vlistDefVarCode(vlistID, varID[i], varCodes[i]);
      // and set the data type
      vlistDefVarDatatype(vlistID, varID[i], datatype);
      // create memory for variables
      var[i] = (double *) Malloc(varSize[i] * sizeof(var[i][0]));
    }

  var_scale(datatype, &mscale, &mrscale);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset
  streamID = composeStream(NULL, prefix, -1, suffix, filetype);

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  {
    uint32_t checksum_state[nvars] = { 0, 0 };
    // Loop over the number of time steps
    for (size_t tsID = 0; tsID < (size_t) nts; tsID++)
      {
        int vdatetime[2] = { 120000, 19850101 + (int) tsID };
        // Set the verification date to 1985-01-01 + tsID
        taxisDefVdate(taxisID, vdatetime[1]);
        // Set the verification time to 12:00:00
        taxisDefVtime(taxisID, vdatetime[0]);
        // Define the time step
        streamDefTimestep(streamID, (int) tsID);

        // Init var1 and var2
        for (size_t j = 0; j < (size_t) nlat; j++)
          for (size_t i = 0; i < (size_t) nlon; i++)
            var[0][i + j * (size_t) nlon] = dg_wobble((double) ((i + tsID) % (size_t) nlon) / (double) (nlon - 1),
                                                      (double) j / (double) (nlat - 1), mscale, mrscale);
        for (size_t k = 0; k < (size_t) nlev; ++k)
          for (size_t j = 0; j < (size_t) nlat; j++)
            for (size_t i = 0; i < (size_t) nlon; i++)
              var[1][i + j * (size_t) nlon + k * (size_t) nlon * (size_t) nlat] = dg_wobble(
                  (double) j / (double) (nlat - 1), (double) ((i + tsID) % (size_t) nlon) / (double) (nlon - 1), mscale, mrscale);

        if (filetype == CDI_FILETYPE_EXT)
          {
            /* EXTRA doesn't store time, only date
             * set the value to 0 before checksumming, because a
             * time field of 0 is what reading an EXTRA file will
             * return */
            vdatetime[0] = 0;
          }
        memcrc_r(&checksum_state[0], (const unsigned char *) vdatetime, sizeof(vdatetime));
        memcrc_r(&checksum_state[0], (const unsigned char *) var[0], varSize[0] * sizeof(var[0][0]));
        memcrc_r(&checksum_state[1], (const unsigned char *) vdatetime, sizeof(vdatetime));
        memcrc_r(&checksum_state[1], (const unsigned char *) var[1], varSize[1] * sizeof(var[1][0]));

        // Write var1 and var2
#ifdef TEST_CHUNK_WRITE
        {
          size_t maxChunkSize = ((size_t) nlon + 1) / 2 * (size_t) nlat * (size_t) nlev;
          double *chunkBuf = (double *) Malloc(maxChunkSize * sizeof(double));
          int varShape[2][3] = { { nlon, nlat, 1 }, { nlon, nlat, nlev } },
              chunk[3][2] = { { 0, nlon / 2 - 1 }, { 0, nlat - 1 }, { 0, 0 } };
          chunk[0][0] = 0;
          chunk[0][1] = nlon / 2 - 1;
          chunk[2][1] = 0;
          get_chunk(chunkBuf, var[0], varShape[0], chunk);
          streamWriteVarChunk(streamID, varID[0], (const int(*)[2]) chunk, chunkBuf, nmiss);
          chunk[2][1] = nlev - 1;
          get_chunk(chunkBuf, var[1], varShape[1], chunk);
          streamWriteVarChunk(streamID, varID[1], (const int(*)[2]) chunk, chunkBuf, nmiss);
          chunk[0][0] = chunk[0][1] + 1;
          chunk[0][1] = nlon - 1;
          chunk[2][1] = 0;
          get_chunk(chunkBuf, var[0], varShape[0], chunk);
          streamWriteVarChunk(streamID, varID[0], (const int(*)[2]) chunk, chunkBuf, nmiss);
          chunk[2][1] = nlev - 1;
          get_chunk(chunkBuf, var[1], varShape[1], chunk);
          streamWriteVarChunk(streamID, varID[1], (const int(*)[2]) chunk, chunkBuf, nmiss);
          free(chunkBuf);
        }
#else
        streamWriteVar(streamID, varID[0], var[0], nmiss);
        streamWriteVar(streamID, varID[1], var[1], nmiss);
#endif
      }
    // write checksums to table file
    {
      FILE *tablefp;
      {
        char *fname = NULL;
        composeFilename(&fname, prefix, -1, "cksum");
        if (!(tablefp = fopen(fname, "w")))
          {
            perror("failed to open table file");
            exit(EXIT_FAILURE);
          }
        free(fname);
      }
      for (size_t i = 0; i < (size_t) nvars; ++i)
        {
          uint32_t cksum;
          int code;
          cksum = memcrc_finish(&checksum_state[i], (off_t) ((varSize[i] * sizeof(var[i][0]) + sizeof(int) * 2) * (size_t) nts));
          code = vlistInqVarCode(vlistID, varID[i]);
          if (fprintf(tablefp, "%08lx %d\n", (unsigned long) cksum, code) < 0)
            {
              perror("failed to write table file");
              exit(EXIT_FAILURE);
            }
        }
      fclose(tablefp);
    }
  }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  for (size_t i = 0; i < nvars; ++i) free(var[i]);
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID[0]);
  zaxisDestroy(zaxisID[1]);
  gridDestroy(gridID);
  return 0;
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
