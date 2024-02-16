#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef USE_MPI
#include <mpi.h>
#include <yaxt.h>
#endif

#include "cdi.h"
#ifdef USE_MPI
#include "cdipio.h"
#include "pio_util.h"
#ifdef HAVE_PPM_CORE
#include <ppm/ppm_uniform_partition.h>
#endif
#endif

#include "cksum.h"
#include "dmemory.h"
#include "error.h"
#include "pio_write.h"
#include "pio_write_setup_grid.h"

#include "simple_model_helper.h"
#include "cdi_uuid.h"

enum
{
  ntfiles = 2,
};

static void
modelRegionCompute(double region[], size_t offset, size_t len, int nlev, int nlat, int nlon, int tsID, double mscale,
                   double mrscale)
{
  size_t local_pos;
  (void) nlev;
  for (local_pos = 0; local_pos < len; ++local_pos)
    {
      size_t global_pos = offset + local_pos;
      int k = (int) (global_pos / (size_t) (nlon * nlat));
      int j = (int) (global_pos % (size_t) (nlon * nlat) / (size_t) nlon);
      int i = (int) (global_pos % (size_t) nlon);
      region[local_pos] = dg_wobble((double) ((i + tsID) % nlon) / (double) (nlon - 1),
                                    (double) ((j + k) % nlat) / (double) (nlat - 1), mscale, mrscale);
    }
}

void
modelRun(const struct model_config *setup, MPI_Comm comm)
{
  struct varDesc_t
  {
    size_t size;
    int nlev, zaxisID, id, code;
    uint32_t checksum_state;
#if USE_MPI
    int chunkSize, start;
    Xt_idxlist partDesc;
#endif
    bool useFloat;
  } * varDesc;
  int gridID, taxisID, vlistID, tsID, tfID = 0;
  enum
  {
    nmiss = 0
  };
  double *levs;
  double *var = NULL, *varslice = NULL;
  float *varsliceF = NULL;
  double mscale, mrscale;
  time_t current_time;
  int vdate = 19850101, vtime = 120000;
  int rank = 0;
  char *filename = NULL;
  int nlon = setup->nlon, nlat = setup->nlat;
  size_t nVars = (size_t) setup->nvars;
  size_t varslice_size = 0, varsliceF_size = 0;
#if USE_MPI
  int *chunks = NULL, *displs = NULL, comm_size = 1;
  Xt_idxlist *partDescPreset = NULL;
  int *conversion = NULL;
#endif

#if USE_MPI
  xmpi(MPI_Comm_rank(comm, &rank));
  xmpi(MPI_Comm_size(comm, &comm_size));
  if (rank == 0 && (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG))
    {
      chunks = Malloc((size_t) comm_size * sizeof(chunks[0]));
      displs = Malloc((size_t) comm_size * sizeof(displs[0]));
      var = Malloc((size_t) nlon * (size_t) nlat * (size_t) setup->max_nlev * sizeof(var[0]));
    }
#else
  (void) comm;
#endif

  var_scale(setup->datatype, &mscale, &mrscale);

  gridID = setupGrid(setup, comm);

  levs = (double *) Malloc((size_t) setup->max_nlev * sizeof(levs[0]));
  {
    double lscale = 1.0 / (double) (setup->max_nlev - 1);
    for (size_t i = 0; i < (size_t) setup->max_nlev; ++i) levs[i] = 101300.0 - 13000.0 * expm1(2.173 * (double) i * lscale);
  }
  vlistID = vlistCreate();

  varDesc = (struct varDesc_t *) Malloc(nVars * sizeof(varDesc[0]));
#if USE_MPI
  if (setup->flags & PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_FLAG)
    {
      partDescPreset = (Xt_idxlist *) Malloc(nVars * sizeof(*partDescPreset));
      conversion = (int *) Malloc(nVars * sizeof(*conversion));
    }
#endif
  for (size_t varIdx = 0; varIdx < nVars; varIdx++)
    {
      int varLevs = (int) cdi_repeatable_random() % 4;
      switch (varLevs)
        {
        case 1: varLevs = setup->max_nlev / 3; break;
        case 2: varLevs = setup->max_nlev >= 11 ? 11 : setup->max_nlev / 2; break;
        case 3: varLevs = setup->max_nlev - 1; break;
        }
      ++varLevs;
      varDesc[varIdx].nlev = varLevs;
      for (size_t i = 0; i < (size_t) varIdx; ++i)
        if (varDesc[i].nlev == varLevs)
          {
            varDesc[varIdx].zaxisID = varDesc[i].zaxisID;
            goto zaxisIDset;
          }
      if (varLevs == 1)
        varDesc[varIdx].zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
      else
        {
          varDesc[varIdx].zaxisID = zaxisCreate(ZAXIS_PRESSURE, varDesc[varIdx].nlev);
          zaxisDefLevels(varDesc[varIdx].zaxisID, levs);
        }
      if (setup->flags & PIO_WRITE_CONFIG_CREATE_UUID_FLAG)
        {
          unsigned char uuid[CDI_UUID_SIZE];
          if (rank == 0) cdiCreateUUID(uuid);
#if USE_MPI
          MPI_Bcast(uuid, CDI_UUID_SIZE, MPI_UNSIGNED_CHAR, 0, comm);
#endif
          cdiDefKeyBytes(varDesc[varIdx].zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);
        }
    zaxisIDset:
      varDesc[varIdx].id = vlistDefVar(vlistID, gridID, varDesc[varIdx].zaxisID, TIME_VARIABLE);
      varDesc[varIdx].size = (size_t) nlon * (size_t) nlat * (size_t) varDesc[varIdx].nlev;
#ifdef USE_MPI
      {
        for (size_t i = 0; i < varIdx; ++i)
          if (varDesc[i].nlev == varLevs)
            {
              varDesc[varIdx].partDesc = varDesc[i].partDesc;
              varDesc[varIdx].start = varDesc[i].start;
              varDesc[varIdx].chunkSize = varDesc[i].chunkSize;
              goto partDescriptionSet;
            }
        struct PPM_extent range = PPM_uniform_partition((struct PPM_extent){ 0, (int32_t) varDesc[varIdx].size }, comm_size, rank);
        int start = range.first;
        int chunkSize = range.size;
        varDesc[varIdx].start = start;
        varDesc[varIdx].chunkSize = chunkSize;
        varDesc[varIdx].partDesc = xt_idxstripes_new(&(struct Xt_stripe){ .start = start, .nstrides = chunkSize, .stride = 1 }, 1);
      }
    partDescriptionSet:;
#endif
      varDesc[varIdx].code = GRIB_USERDEF + (int) varIdx;
      vlistDefVarCode(vlistID, varDesc[varIdx].id, varDesc[varIdx].code);
      vlistDefVarDatatype(vlistID, varDesc[varIdx].id, setup->datatype);
      varDesc[varIdx].useFloat = (bool) (cdi_repeatable_random() & 1);
#ifdef USE_MPI
      if (setup->flags & PIO_WRITE_CONFIG_PRESET_DECOMPOSITION_FLAG)
        {
          partDescPreset[varIdx] = varDesc[varIdx].partDesc;
          conversion[varIdx] = varDesc[varIdx].useFloat ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
        }
#endif
    }

  taxisID = taxisCreate(setup->taxistype);
  if (setup->taxisunit != -1) taxisDefTunit(taxisID, setup->taxisunit);
  vlistDefTaxis(vlistID, taxisID);

  for (tfID = 0; tfID < ntfiles; tfID++)
    {
      for (size_t varIdx = 0; varIdx < nVars; ++varIdx) varDesc[varIdx].checksum_state = 0;

      int streamID = composeStream(&filename, setup->prefix, tfID, setup->suffix, setup->filetype);
#ifdef USE_MPI
      if (partDescPreset)
        cdiPioStreamDefDecomposedVlist(streamID, vlistID, partDescPreset, conversion);
      else
#endif
        streamDefVlist(streamID, vlistID);

      vdate = 19850101;
      vtime = 120000;
      current_time = cditime2time_t(vdate, vtime);
      for (tsID = 0; tsID < setup->nts; tsID++)
        {
          int vdatetime[2];
          time_t2cditime(current_time, &vdatetime[1], &vdatetime[0]);
          taxisDefVdate(taxisID, vdatetime[1]);
          taxisDefVtime(taxisID, vdatetime[0]);
          streamDefTimestep(streamID, tsID);
          if (setup->filetype == CDI_FILETYPE_EXT)
            {
              /* EXTRA doesn't store time, only date
               * set the value to 0 before checksumming, because a
               * time field of 0 is what reading an EXTRA file will
               * return */
              vdatetime[0] = 0;
            }
          for (size_t varIdx = 0; varIdx < nVars; ++varIdx)
            {
#ifdef USE_MPI
              int start = varDesc[varIdx].start;
              size_t chunkSize = (size_t) varDesc[varIdx].chunkSize;
#else
              size_t chunkSize = varDesc[varIdx].size;
              int start = 0;
#endif
              if (varslice_size < chunkSize)
                {
                  varslice = (double *) Realloc(varslice, chunkSize * sizeof(var[0]));
                  varslice_size = chunkSize;
                }
              modelRegionCompute(varslice, (size_t) start, chunkSize, varDesc[varIdx].nlev, nlat, nlon, tsID, mscale, mrscale);
              bool useFloat = varDesc[varIdx].useFloat;
              if (useFloat)
                {
                  if (varsliceF_size < chunkSize)
                    {
                      varsliceF = (float *) Realloc(varsliceF, chunkSize * sizeof(varsliceF[0]));
                      varsliceF_size = chunkSize;
                    }
                  for (size_t i = 0; i < chunkSize; ++i) varslice[i] = varsliceF[i] = (float) varslice[i];
                }

              if (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG)
                {
#if USE_MPI
                  int chunk = (int) chunkSize;
                  xmpi(MPI_Gather(&chunk, 1, MPI_INT, chunks, 1, MPI_INT, 0, comm));
                  if (rank == 0)
                    {
                      int accum = 0;
                      for (size_t i = 0; i < (size_t) comm_size; ++i)
                        {
                          displs[i] = accum;
                          accum += chunks[i];
                        }
                    }
                  xmpi(MPI_Gatherv(varslice, chunk, MPI_DOUBLE, var, chunks, displs, MPI_DOUBLE, 0, comm));
#else
                  var = varslice;
#endif
                }
              if (rank == 0 && (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG))
                {
                  memcrc_r(&varDesc[varIdx].checksum_state, (const unsigned char *) vdatetime, sizeof(vdatetime));
                  memcrc_r(&varDesc[varIdx].checksum_state, (const unsigned char *) var, varDesc[varIdx].size * sizeof(var[0]));
                }

#ifdef USE_MPI
              if (useFloat)
                streamWriteVarPartF(streamID, varDesc[varIdx].id, varsliceF, nmiss, varDesc[varIdx].partDesc);
              else
                streamWriteVarPart(streamID, varDesc[varIdx].id, varslice, nmiss, varDesc[varIdx].partDesc);
#else
              if (useFloat)
                streamWriteVarF(streamID, varDesc[varIdx].id, varsliceF, nmiss);
              else
                streamWriteVar(streamID, varDesc[varIdx].id, varslice, nmiss);
#endif
            }
          current_time += 86400;
#ifdef USE_MPI
          pioWriteTimestep();
#endif
        }
      streamClose(streamID);
      if (rank == 0 && (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG))
        {
          FILE *tablefp;
          composeFilename(&filename, setup->prefix, tfID, "cksum");
          if (!(tablefp = fopen(filename, "w")))
            {
              perror("failed to open table file");
              exit(EXIT_FAILURE);
            }
          for (size_t varIdx = 0; varIdx < nVars; ++varIdx)
            {
              uint32_t cksum
                  = memcrc_finish(&varDesc[varIdx].checksum_state,
                                  (off_t) ((varDesc[varIdx].size * sizeof(var[0]) + sizeof(int) * 2) * (size_t) setup->nts));
              int code = vlistInqVarCode(vlistID, varDesc[varIdx].id);
              if (fprintf(tablefp, "%08lx %d\n", (unsigned long) cksum, code) < 0)
                {
                  perror("failed to write table file");
                  exit(EXIT_FAILURE);
                }
            }
          if (fclose(tablefp) != 0)
            {
              perror("failed to close table file");
              exit(EXIT_FAILURE);
            }
        }
    }
#ifdef USE_MPI
  pioEndTimestepping();
#endif
  Free(varslice);
  Free(varsliceF);
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  for (size_t varIdx = 0; varIdx < nVars; varIdx++)
    {
      int zID = varDesc[varIdx].zaxisID;
      if (zID != CDI_UNDEFID)
        {
          zaxisDestroy(zID);
#if USE_MPI
          xt_idxlist_delete(varDesc[varIdx].partDesc);
#endif
          for (size_t j = varIdx + 1; j < nVars; ++j)
            if (zID == varDesc[j].zaxisID) varDesc[j].zaxisID = CDI_UNDEFID;
        }
    }
  gridDestroy(gridID);
#if USE_MPI
  Free(conversion);
  Free(partDescPreset);
  Free(displs);
  Free(chunks);
  Free(var);
#endif
  Free(varDesc);
  Free(levs);
  Free(filename);
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
