#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
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
#include <core/ppm_combinatorics.h>
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
  nproma = 16,
};

static void
modelRegionCompute(double region[], int nlev, int nlat, int nlon, const int chunkStart[3], const int chunkSize[3], int tsID,
                   double mscale, double mrscale)
{
  (void) nlev;
  int is = chunkStart[0], js = chunkStart[1], ks = chunkStart[2], m = chunkSize[0], n = chunkSize[1], o = chunkSize[2],
      jstride = chunkSize[0], kstride = ((jstride * chunkSize[1] + nproma - 1) / nproma) * nproma;

  tsID %= nlon;
  for (int k = 0; k < o; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i)
        region[k * kstride + j * jstride + i]
            = dg_wobble((double) ((i + is + tsID) % nlon) / (double) (nlon - 1),
                        (double) ((j + js + k + ks) % nlat) / (double) (nlat - 1), mscale, mrscale);
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
    int chunkSize[2], start[2];
    Xt_idxlist partDesc;
    Xt_redist redist4gather;
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
#if !USE_MPI
  float *varF = NULL;
#endif
  double mscale, mrscale;
  time_t current_time;
  int vdate = 19850101, vtime = 120000;
  int rank = 0;
  char *filename = NULL;
  int nlon = setup->nlon, nlat = setup->nlat;
  size_t nVars = setup->nvars > 0 ? (size_t) setup->nvars : 0;
  size_t varslice_size = 0, varsliceF_size = 0;
  ;
#if USE_MPI
  int comm_size = 1;
  int npart[2], rank_coord[2];
  int *blk_displ, *blk_lens;
#endif

#if USE_MPI
  xmpi(MPI_Comm_rank(comm, &rank));
  xmpi(MPI_Comm_size(comm, &comm_size));
#else
  (void) comm;
#endif

#if USE_MPI
  bool needsGather = setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG;
#else
  bool needsGather = true;
#endif
  if (rank == 0 && (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG))
    {
      var = (double *) Malloc((size_t) nlon * (size_t) nlat * (size_t) setup->max_nlev * sizeof(var[0]));
#if !USE_MPI
      varF = (float *) Malloc((size_t) nlon * (size_t) nlat * (size_t) setup->max_nlev * sizeof(varF[0]));
#endif
    }

#if USE_MPI
  if (comm_size == 1)
    {
      npart[0] = 1;
      npart[1] = 1;
      rank_coord[0] = 0;
      rank_coord[1] = 0;
    }
  else
    {
      findPartition2D(npart, comm_size);
      rank_coord[0] = rank % npart[0], rank_coord[1] = rank / npart[0];
    }
  blk_displ = Malloc((size_t) setup->max_nlev * sizeof(blk_displ[0]) * 2);
  blk_lens = blk_displ + setup->max_nlev;
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
              varDesc[varIdx].redist4gather = varDesc[i].redist4gather;
              varDesc[varIdx].partDesc = varDesc[i].partDesc;
              for (size_t j = 0; j < 2; ++j)
                {
                  varDesc[varIdx].start[j] = varDesc[i].start[j];
                  varDesc[varIdx].chunkSize[j] = varDesc[i].chunkSize[j];
                }
              goto partDescriptionSet;
            }
        int start[2], chunkSize[3], varSize[2] = { nlon, nlat };
        for (size_t i = 0; i < 2; ++i)
          {
            struct PPM_extent range = PPM_uniform_partition((struct PPM_extent){ 0, varSize[i] }, npart[i], rank_coord[i]);
            start[i] = range.first;
            chunkSize[i] = range.size;
            fprintf(stderr, "%d: start[%zu]=%d, chunkSize[%zu] = %d\n", rank, i, start[i], i, chunkSize[i]);
            varDesc[varIdx].start[i] = range.first;
            varDesc[varIdx].chunkSize[i] = range.size;
          }
        Xt_int varSizeXt[3] = { (Xt_int) nlon, (Xt_int) nlat, (Xt_int) varLevs };
        chunkSize[2] = varLevs;
        Xt_int varStartXt[3] = { start[0], start[1], 0 };
        Xt_idxlist part_idxlist = xt_idxsection_new(0, (varLevs > 1 ? 3 : 2), varSizeXt, chunkSize, varStartXt), gather_idxlist;
        varDesc[varIdx].partDesc = part_idxlist;
        if (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG)
          {
            if (rank == 0)
              {
                gather_idxlist
                    = xt_idxstripes_new(&(struct Xt_stripe){ .start = 0, .stride = 1, .nstrides = (int) varDesc[varIdx].size }, 1);
              }
            else
              gather_idxlist = xt_idxempty_new();
            Xt_xmap xmap4gather = xt_xmap_all2all_new(part_idxlist, gather_idxlist, comm);
            xt_idxlist_delete(gather_idxlist);
            struct Xt_offset_ext *src_blocks = Malloc((size_t) varLevs * sizeof(*src_blocks));
            struct Xt_offset_ext dst_block = { .start = 0, .size = nlon * nlat * varLevs, .stride = 1 };
            size_t levStride = (((size_t) chunkSize[0] * (size_t) chunkSize[1] + nproma - 1) / nproma) * nproma;
            for (size_t i = 0; i < (size_t) varLevs; ++i)
              src_blocks[i]
                  = (struct Xt_offset_ext){ .start = (int) (i * levStride), .size = chunkSize[0] * chunkSize[1], .stride = 1 };
            varDesc[varIdx].redist4gather = xt_redist_p2p_ext_new(xmap4gather, varLevs, src_blocks, 1, &dst_block, MPI_DOUBLE);
            Free(src_blocks);
            xt_xmap_delete(xmap4gather);
          }
      }
    partDescriptionSet:;
#endif
      varDesc[varIdx].code = GRIB_USERDEF + (int) varIdx;
      vlistDefVarCode(vlistID, varDesc[varIdx].id, varDesc[varIdx].code);
      vlistDefVarDatatype(vlistID, varDesc[varIdx].id, setup->datatype);
      varDesc[varIdx].useFloat = (bool) (cdi_repeatable_random() & 1);
    }

  taxisID = taxisCreate(setup->taxistype);
  vlistDefTaxis(vlistID, taxisID);

  for (tfID = 0; tfID < ntfiles; tfID++)
    {
      for (size_t varIdx = 0; varIdx < nVars; ++varIdx) varDesc[varIdx].checksum_state = 0;

      int streamID = composeStream(&filename, setup->prefix, tfID, setup->suffix, setup->filetype);
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
              size_t varLevs = (size_t) varDesc[varIdx].nlev;
#ifdef USE_MPI
              int start[3] = { varDesc[varIdx].start[0], varDesc[varIdx].start[1], 0 };
              int chunk[3] = { varDesc[varIdx].chunkSize[0], varDesc[varIdx].chunkSize[1], (int) varLevs };
#else
              int chunk[3] = { nlon, nlat, (int) varLevs };
              int start[3] = { 0, 0, 0 };
#endif
              size_t chunkSize
                  = (((size_t) chunk[0] * (size_t) chunk[1] + (size_t) (nproma - 1)) / (size_t) nproma) * (size_t) nproma * varLevs;
              if (varslice_size < chunkSize)
                {
                  varslice = (double *) Realloc(varslice, chunkSize * sizeof(var[0]));
                  varslice_size = chunkSize;
                }
              modelRegionCompute(varslice, (int) varLevs, nlat, nlon, start, chunk, tsID, mscale, mrscale);
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

              if (needsGather)
                {
#if USE_MPI
                  xt_redist_s_exchange1(varDesc[varIdx].redist4gather, varslice, var);
#else
                  size_t layerSize = (size_t) (chunk[0] * chunk[1]);
                  size_t nblk = (layerSize + nproma - 1) / nproma - 1;
                  size_t npromz = layerSize - nblk * nproma;
                  for (size_t k = 0; k < varLevs; ++k)
                    {
                      for (size_t j = 0; j < nblk; ++j)
                        for (size_t i = 0; i < nproma; ++i)
                          var[k * layerSize + j * nproma + i] = varslice[k * (nblk + 1) * nproma + j * nproma + i];
                      for (size_t i = 0; i < npromz; ++i)
                        var[k * layerSize + nblk * nproma + i] = varslice[k * (nblk + 1) * nproma + nblk * nproma + i];
                    }
                  if (useFloat)
                    for (size_t k = 0; k < varLevs; ++k)
                      for (size_t i = 0; i < layerSize; ++i) varF[k * layerSize + i] = (float) var[k * layerSize + i];
#endif
                }
              if (rank == 0 && (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG))
                {
                  memcrc_r(&varDesc[varIdx].checksum_state, (const unsigned char *) vdatetime, sizeof(vdatetime));
                  memcrc_r(&varDesc[varIdx].checksum_state, (const unsigned char *) var, varDesc[varIdx].size * sizeof(var[0]));
                }

#ifdef USE_MPI
              size_t layerSize = (size_t) (chunk[0] * chunk[1]);
              size_t nblk = (layerSize + nproma - 1) / nproma - 1;
              for (size_t k = 0; k < varLevs; ++k)
                {
                  blk_displ[k] = (int) (k * (nblk + 1) * nproma);
                  blk_lens[k] = (int) layerSize;
                }
              if (useFloat)
                streamWriteScatteredVarPartF(streamID, varDesc[varIdx].id, varsliceF, (int) varLevs, blk_lens, blk_displ, nmiss,
                                             varDesc[varIdx].partDesc);
              else
                streamWriteScatteredVarPart(streamID, varDesc[varIdx].id, varslice, (int) varLevs, blk_lens, blk_displ, nmiss,
                                            varDesc[varIdx].partDesc);
#else
              if (useFloat)
                streamWriteVarF(streamID, varDesc[varIdx].id, varF, nmiss);
              else
                streamWriteVar(streamID, varDesc[varIdx].id, var, nmiss);
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
              uint32_t cksum;
              int code;
              cksum = memcrc_finish(&varDesc[varIdx].checksum_state,
                                    (off_t) ((varDesc[varIdx].size * sizeof(var[0]) + sizeof(int) * 2) * (size_t) setup->nts));
              code = vlistInqVarCode(vlistID, varDesc[varIdx].id);
              if (fprintf(tablefp, "%08lx %d\n", (unsigned long) cksum, code) < 0)
                {
                  perror("failed to write table file");
                  exit(EXIT_FAILURE);
                }
            }
          fclose(tablefp);
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
          if (setup->flags & PIO_WRITE_CONFIG_CHECKSUM_FLAG) xt_redist_delete(varDesc[varIdx].redist4gather);
#endif
          for (size_t j = varIdx + 1; j < nVars; ++j)
            if (zID == varDesc[j].zaxisID) varDesc[j].zaxisID = CDI_UNDEFID;
        }
    }
  gridDestroy(gridID);
  Free(var);
#if USE_MPI
  Free(blk_displ);
#else
  Free(varF);
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
