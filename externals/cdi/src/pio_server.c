/** @file ioServer.c
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#ifdef HAVE_PARALLEL_NC4
#include <setjmp.h>
#endif
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_PARALLEL_NC4
#include <core/ppm_extents.h>
#include <core/ppm_combinatorics.h>
#include <core/ppm_rectilinear.h>
#include <ppm/ppm_uniform_partition.h>
#endif
#include <mpi.h>
#include <yaxt.h>

#include "cdi.h"
#include "cdipio.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "namespace.h"
#include "taxis.h"
#include "pio.h"
#include "pio_cdf_int.h"
#include "pio_comm.h"
#include "pio_conf.h"
#include "pio_idxlist_cache.h"
#include "pio_id_set.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_server.h"
#include "pio_util.h"
#include "pio_xmap_cache.h"
#ifndef HAVE_NETCDF_PAR_H
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-macros"
#endif
/* the following causes netcdf.h to not try to redefine MPI_Comm etc. */
#define MPI_INCLUDED
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
#endif
#include "resource_handle.h"
#include "resource_unpack.h"
#include "stream_cdf.h"
#include "vlist_var.h"

struct clientBuf
{
  size_t size;
  unsigned char *mem;
  int dictSize;
};

struct streamMemLayout
{
  Xt_uid varPartIdxListUID;
  int offset, conversion;
};

struct cacheRedist
{
  Xt_redist redist;
  size_t sliceSize;
#if defined HAVE_PARALLEL_NC4
  // additional member for parallel writing variant for NetCDF4
  struct PPM_extent sliceExtent[3];
#endif
};

static void *sharedClientBuf = NULL;
static size_t sharedClientBufSize = 0;

static struct
{
  MPI_Win getWin;
  struct clientBuf *clientBuf;
#ifdef HAVE_LIBNETCDF
  int ownerRank;
#endif
  /* put data for description of last layout from RMA GET here */
  struct streamMemLayout *prevLayout;
  size_t numRetained;
  struct cacheRedist *retained;
  struct partDescPreset clientDeco;
} * rxWin;

static struct
{
  size_t size, used;
  void *mem;
} aggBuf;

/* map streamID's to linear index into arrays */
static struct idList openStreams;

static struct cdiPioIdxlistCache *DstIdxlistCache;
static size_t neededDstIdxlistCacheSize;

static struct xmapCache *XmapCache;
static size_t neededXmapCacheSize;

struct recordWrite
{
  int varID, level;
  size_t dataSize;
};

struct streamMapping
{
  /* data entry varMap[i] contains data for variable i or -1 if no
   * data entry for i has been transferred */
  int *varMap;
  /* numLvls[i] number of levels written for variable i or 0 if
   * variable is not written to this timestep */
  int *numLvlsW;
  /* nMiss[i] = missing values were provided for variable i */
  int *hasMissing;
  int numWrittenRecords;
  int numVars;
  struct streamMemLayout *layout;
  struct recordWrite writtenRecords[];
};

#ifdef HAVE_PARALLEL_NC4
/* prime factorization of number of pio collectors */
static uint32_t *pioPrimes;
static int numPioPrimes;
#endif

/************************************************************************/

static void
cdiPioServerStreamWinDestroy(size_t streamIdx, const struct cdiPioConf *conf)
{
  if (rxWin[streamIdx].getWin != MPI_WIN_NULL)
    {
      if (conf->batchedRMA) Free(rxWin[streamIdx].clientBuf[0].mem);
      rxWin[streamIdx].clientBuf[0].mem = NULL;
      xmpi(MPI_Win_free(&rxWin[streamIdx].getWin));
    }
}

static int numClients_, *clientRanks_;

static void
setupClientRanks(void)
{
  MPI_Group clientGroup = cdiPioInqRemoteGroup();
  xmpi(MPI_Group_size(clientGroup, &numClients_));
  clientRanks_ = Malloc((size_t) numClients_ * sizeof(clientRanks_[0]));
  int *ranks = Malloc((size_t) numClients_ * sizeof(ranks[0]));
  for (int i = 0; i < numClients_; ++i) ranks[i] = i;
  MPI_Comm collClientIntraComm = cdiPioInqCollClientIntraComm();
  MPI_Group groupCollClient;
  xmpi(MPI_Comm_group(collClientIntraComm, &groupCollClient));
  xmpi(MPI_Group_translate_ranks(clientGroup, numClients_, ranks, groupCollClient, clientRanks_));
  xmpi(MPI_Group_free(&groupCollClient));
  Free(ranks);
}

static void
createClientStreamBuf(size_t streamIdx, const struct clientBufSize *bufSizes, const struct cdiPioConf *conf)
{
  /* find and tabulate aggregate size needed for all clients of collector */
  size_t streamBufferSize = 0;
  for (size_t i = 0; i < (size_t) numClients_; ++i)
    {
      streamBufferSize += (rxWin[streamIdx].clientBuf[i].size = bufSizes[i].bufSize);
      rxWin[streamIdx].clientBuf[i].dictSize = bufSizes[i].numDataRecords + bufSizes[i].numRPCRecords;
    }
  /* set pointer to RMA buffer for client 0 of collector */
  if (conf->batchedRMA)
    rxWin[streamIdx].clientBuf[0].mem = Malloc(streamBufferSize);
  else
    {
      if (streamBufferSize > sharedClientBufSize)
        {
          intptr_t oldmem = (intptr_t) sharedClientBuf;
          sharedClientBuf = Realloc(sharedClientBuf, streamBufferSize);
          sharedClientBufSize = streamBufferSize;
          if (oldmem != (intptr_t) sharedClientBuf)
            {
              size_t numStreams = openStreams.size;
              for (size_t j = 0; j < numStreams; ++j)
                if (rxWin[j].getWin != MPI_WIN_NULL)
                  {
                    unsigned char *newmem = sharedClientBuf;
                    for (size_t i = 0; i < (size_t) numClients_; ++i)
                      {
                        rxWin[j].clientBuf[i].mem = newmem;
                        newmem += rxWin[j].clientBuf[i].size;
                      }
                  }
            }
        }
      rxWin[streamIdx].clientBuf[0].mem = sharedClientBuf;
    }
  /* set pointers for other clients */
  {
    unsigned char *newmem = rxWin[streamIdx].clientBuf[0].mem;
    for (size_t i = 1; i < (size_t) numClients_; ++i)
      {
        newmem += rxWin[streamIdx].clientBuf[i - 1].size;
        rxWin[streamIdx].clientBuf[i].mem = newmem;
      }
  }
}

static void
cdiPioServerStreamWinCreate(size_t streamIdx, MPI_Comm collClientIntraComm)
{
  MPI_Info no_locks_info;
  xmpi(MPI_Info_create(&no_locks_info));
  xmpi(MPI_Info_set(no_locks_info, "no_locks", "true"));
  xmpi(MPI_Win_create(MPI_BOTTOM, 0, 1, no_locks_info, collClientIntraComm, &rxWin[streamIdx].getWin));
  xmpi(MPI_Info_free(&no_locks_info));
}

/************************************************************************/

static void
readFuncCall(struct winHeaderEntry *header, size_t streamIdx)
{
  int funcID = header->id;
  union funcArgs *funcArgs = &(header->specific.funcArgs);

  xassert(funcID >= MINFUNCID && funcID <= MAXFUNCID);
  switch (funcID)
    {
    case STREAMDEFTIMESTEP:
      {
        MPI_Comm pioInterComm = cdiPioInqInterComm();
        int streamID = funcArgs->streamNewTimestep.streamID;
        int originNamespace = namespaceResHDecode(streamID).nsp;
        streamID = namespaceAdaptKey2(streamID);
        int oldTaxisID = vlistInqTaxis(streamInqVlist(streamID));
        int position = header->offset;
        int changedTaxisID = taxisUnpack((char *) rxWin[streamIdx].clientBuf[0].mem, (int) rxWin[streamIdx].clientBuf[0].size,
                                         &position, originNamespace, &pioInterComm, 0);
        taxis_t *oldTaxisPtr = taxisPtr(oldTaxisID);
        taxis_t *changedTaxisPtr = taxisPtr(changedTaxisID);
        ptaxisCopy(oldTaxisPtr, changedTaxisPtr);
        taxisDestroy(changedTaxisID);
        streamDefTimestep(streamID, funcArgs->streamNewTimestep.tsID);
      }
      break;
    default: xabort("REMOTE FUNCTIONCALL NOT IMPLEMENTED!");
    }
}

/************************************************************************/

static void
resizeVarGatherBuf(size_t size, void **buf, size_t *bufSize)
{
  if (size <= *bufSize)
    ;
  else
    *buf = Realloc(*buf, *bufSize = size);
}

#define wHECast(buf) ((struct winHeaderEntry *) (void *) buf)

static Xt_xmap
buildVarXmap(struct Xt_offset_ext *restrict partExts, const struct clientBuf *restrict clientBuf, size_t headerIdx,
             Xt_idxlist dstList, Xt_idxlist *partDescPreset, MPI_Comm pioInterComm, MPI_Comm collComm, int varID,
             const struct cdiPioConf *conf)
{
  size_t numClients = (size_t) numClients_;
  Xt_idxlist *part = partDescPreset ? partDescPreset : Malloc(numClients * sizeof(part[0]));
  int conversion = (wHECast(clientBuf[0].mem))[headerIdx].id;
  size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
  for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx)
    {
      unsigned char *clientMem = clientBuf[clientIdx].mem;
      struct dataRecord *dataHeader = &(wHECast(clientMem))[headerIdx].specific.dataRecord;
      xassert(dataHeader->varID == varID && (wHECast(clientMem))[headerIdx].id == conversion
              && ((wHECast(clientMem))[headerIdx + 1].id == PARTDESCMARKER));
      if (!partDescPreset)
        {
          int position = (wHECast(clientMem))[headerIdx + 1].offset;
          xassert(position > 0 && ((size_t) position >= sizeof(struct winHeaderEntry) * (size_t) clientBuf[clientIdx].dictSize)
                  && ((size_t) position < clientBuf[clientIdx].size));
          part[clientIdx] = xt_idxlist_unpack(clientMem, (int) clientBuf[clientIdx].size, &position, pioInterComm);
        }
      unsigned partSize = (unsigned) xt_idxlist_get_num_indices(part[clientIdx]);
      size_t charOfs = (size_t) ((clientMem + (wHECast(clientMem))[headerIdx].offset) - clientBuf[0].mem);
      xassert(charOfs % elemSize == 0 && charOfs / elemSize + partSize <= INT_MAX);
      int elemOfs = (int) (charOfs / elemSize);
      partExts[clientIdx].start = elemOfs;
      partExts[clientIdx].size = (int) partSize;
      partExts[clientIdx].stride = 1;
    }
  Xt_idxlist srcList = xt_idxlist_collection_new(part, (int) numClients);
  if (!partDescPreset)
    {
      for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx) xt_idxlist_delete(part[clientIdx]);
      Free(part);
    }
  if (conf->stripify)
    {
      Xt_idxlist srcListStriped = xt_idxstripes_from_idxlist_new(srcList);
      xt_idxlist_delete(srcList);
      srcList = srcListStriped;
    }
  Xt_xmap gatherXmap = conf->xmap_new(srcList, dstList, collComm);
  xt_idxlist_delete(srcList);
  return gatherXmap;
}

#ifdef HAVE_PARALLEL_NC4
/*
 * Given a divisor @a composite_div,
 * find subset $F$ of prime divisors of @a composite_div such that the
 * product of $F$ is less than npartMax and maximal
 */
unsigned
findMaxDivision(unsigned npartMax, unsigned composite_div)
{
  uint32_t factors[31], *factors_ = factors;
  int numFactors = PPM_prime_factorization_32((uint32_t) composite_div, &factors_);
  /* try to use prime factors */
  uint_fast32_t divAttempt, maxDiv = 1;
  /* test all possible assignments of factors, starting with
   * only one assigned (omitting no assigned case because that would
   * never be better than start value of maxDiv */
  for (int numAssigned = 1; numAssigned <= numFactors; ++numAssigned)
    {
      uint_fast32_t pattern = (UINT32_C(1) << numAssigned) - 1, lastPattern = pattern << (numFactors - numAssigned);
      do
        {
          divAttempt = 1;
          /* loop over all factors */
          for (uint_fast32_t i = 0; i < (uint_fast32_t) numFactors; ++i)
            {
              uint_fast32_t assigned = (pattern >> i) & 1;
              if (assigned) divAttempt *= factors[i];
            }
          if (divAttempt <= npartMax && divAttempt > maxDiv) maxDiv = divAttempt;
          /* find next sequence of numAssigned set bits and numFactors
           * - numFactors unset bits */
          {
            uint_fast32_t t;
#if HAVE_DECL___BUILTIN_CTZ
            t = pattern | (pattern - 1);
            pattern = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz((unsigned) pattern) + 1));
#else
            t = (pattern | (pattern - 1)) + 1;
            pattern = t | ((((t & -t) / (pattern & -pattern)) >> 1) - 1);
#endif
          }
        }
      while (pattern <= lastPattern);
    }
  return (unsigned) maxDiv;
}

static void
queryVarBounds(struct PPM_extent varShape[3], int vlistID, int varID)
{
  int sizes[3];
  cdiPioQueryVarDims(sizes, vlistID, varID);
  for (unsigned i = 0; i < 3; ++i) varShape[i].first = 0, varShape[i].size = sizes[i];
}

struct xyzDims
{
  unsigned sizes[3];
};

static struct xyzDims
varDimsCollGridMax(const struct PPM_extent varDims[3])
{
  struct xyzDims collGrid = { { 1U, 1U, 1U } };
  unsigned collDiv = (unsigned) (commInqSizeColl());
  for (size_t i = 3; i > 0; --i)
    {
      unsigned usedDiv = collGrid.sizes[i - 1] = findMaxDivision((unsigned) varDims[i - 1].size, collDiv);
      collDiv /= usedDiv;
    }
  return collGrid;
}

/* compute distribution of collectors such that number of collectors
 * <= number of variable grid cells in each dimension */
static struct xyzDims
varDimsCollGridMatch(const struct PPM_extent varDims[3])
{
  struct xyzDims collGrid = { { 1, 1, 1 } };
  /* because of storage order, dividing dimension 3 first is preferred */
  for (int i = 0; i < numPioPrimes; ++i)
    {
      for (int dim = 2; dim >= 0; --dim)
        if (collGrid.sizes[dim] * pioPrimes[i] <= (unsigned) varDims[dim].size)
          {
            collGrid.sizes[dim] *= pioPrimes[i];
            goto nextPrime;
          }
      /* no easy I/O decomposition found, do exhaustive search for not
       * necessarily perfect decomposition */
      return varDimsCollGridMax(varDims);
    nextPrime:;
    }
  return collGrid;
}

static void
myVarPart(struct PPM_extent varShape[3], struct xyzDims collGrid, struct PPM_extent myPart[3])
{
  int32_t myCollGridCoord[3];
  struct PPM_extent collGridShape[3];
  unsigned collGridSize = 1;
  for (size_t i = 0; i < 3; ++i)
    {
      collGridShape[i].first = 0;
      collGridShape[i].size = (int32_t) collGrid.sizes[i];
      collGridSize *= collGrid.sizes[i];
    }
  int collRank = commInqRankColl();
  if ((unsigned) collRank < collGridSize)
    {
      PPM_lidx2rlcoord_e(3, collGridShape, collRank, myCollGridCoord);
      xdebug("my coord: (%d, %d, %d)", myCollGridCoord[0], myCollGridCoord[1], myCollGridCoord[2]);
      int32_t sizes[3];
      for (size_t i = 0; i < 3; ++i) sizes[i] = (int32_t) collGrid.sizes[i];
      PPM_uniform_partition_nd(3, varShape, sizes, myCollGridCoord, myPart);
    }
  else
    for (size_t i = 0; i < 3; ++i) myPart[i].first = myPart[i].size = 0;
}

static void
cdiPioNetCDFParChunk(int vlistID, int varID, Xt_idxlist *preWriteChunk, struct PPM_extent varChunk[3])
{
  struct PPM_extent varShape[3];
  queryVarBounds(varShape, vlistID, varID);
  const struct xyzDims collGrid = varDimsCollGridMatch(varShape);
  xdebug("writing varID %d with dimensions: "
         "x=%d, y=%d, z=%d,\n"
         "found distribution with dimensions:"
         " x=%d, y=%d, z=%d.",
         varID, varShape[0].size, varShape[1].size, varShape[2].size, collGrid.sizes[0], collGrid.sizes[1], collGrid.sizes[2]);
  myVarPart(varShape, collGrid, varChunk);
  {
    Xt_int preWriteChunkStart[3];
    int preWriteChunkSize[3];
    Xt_int varDims[3];
    for (int i = 0; i < 3; ++i)
      {
        varDims[2 - i] = varShape[i].size;
        preWriteChunkStart[2 - i] = (Xt_int) varChunk[i].first;
        preWriteChunkSize[2 - i] = (int) varChunk[i].size;
      }
    *preWriteChunk = cdiPioIdxlistCacheAddSection3D(DstIdxlistCache, varDims, preWriteChunkStart, preWriteChunkSize);
  }
}
#endif

static void
allocUIDLookup(size_t numClients, Xt_uid *restrict *uids, int *restrict *partSizes)
{
  size_t uidBytes = sizeof(**uids) * (numClients + 1), partSizesBytes = sizeof(**partSizes) * (numClients + 1),
         partSizeAlign = sizeof(**partSizes), uidBytesRoundUp = ((uidBytes + partSizeAlign - 1) / partSizeAlign * partSizeAlign);
  *uids = Malloc(uidBytes + uidBytesRoundUp + partSizesBytes);
  *partSizes = (int *) (void *) ((unsigned char *) *uids + uidBytes + uidBytesRoundUp);
}

static Xt_idxlist buildVarSlicesIdxList(int vlistID, int varID, int startLvl, int numLvl, const struct cdiPioConf *conf);

static void
buildDecoPresetXmapsCommon(Xt_idxlist outPartList, size_t numClients, const Xt_uid *partDescUIDs, Xt_idxlist *partDesc,
                           int *restrict partSizes, Xt_uid *uids, MPI_Comm collComm, const struct cdiPioConf *conf)
{
  uids[0] = xt_idxlist_get_uid(outPartList);
  for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx) uids[clientIdx + 1] = partDescUIDs[clientIdx];
  Xt_xmap varXmap;
  if (!(varXmap = cdiPioXmapCacheLookup(XmapCache, uids, partSizes)))
    {
      Xt_idxlist srcList = xt_idxlist_collection_new(partDesc, (int) numClients);
      if (conf->stripify)
        {
          Xt_idxlist srcListStriped = xt_idxstripes_from_idxlist_new(srcList);
          xt_idxlist_delete(srcList);
          srcList = srcListStriped;
        }
      varXmap = conf->xmap_new(srcList, outPartList, collComm);
      xt_idxlist_delete(srcList);
      partSizes[0] = xt_idxlist_get_num_indices(outPartList);
      for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx)
        partSizes[clientIdx + 1] = xt_idxlist_get_num_indices(partDesc[clientIdx]);
      cdiPioXmapCacheAdd(XmapCache, uids, partSizes, varXmap);
    }
}

#ifdef HAVE_LIBNETCDF
typedef Xt_idxlist (*outDstListConstruct)(int vlistID, int varID, const struct cdiPioConf *conf);
#ifdef HAVE_PARALLEL_NC4
Xt_idxlist
createOutDstListNetCDFPar(int vlistID, int varID, const struct cdiPioConf *conf)
{
  (void) conf;
  struct PPM_extent varChunk[3];
  Xt_idxlist preWriteChunk;
  cdiPioNetCDFParChunk(vlistID, varID, &preWriteChunk, varChunk);
  return preWriteChunk;
}
#endif

static Xt_idxlist
createOutDstListNetCDFSerialFunnel(int vlistID, int varID, const struct cdiPioConf *conf)
{
  Xt_idxlist preWriteChunk = buildVarSlicesIdxList(vlistID, varID, -1, -1, conf);
  return preWriteChunk;
}

static Xt_idxlist
createOutDstListNetCDFSerialAssist(int vlistID, int varID, const struct cdiPioConf *conf)
{
  (void) vlistID;
  (void) varID;
  (void) conf;
  return xt_idxempty_new();
}

static void
buildDecoPresetXmapsNetCDF(int vlistID, struct partDescPreset clientDeco, MPI_Comm collComm, outDstListConstruct getOutDstList,
                           const struct cdiPioConf *conf)
{
  int nVars = vlistNvars(vlistID);
  size_t numClients = (size_t) numClients_;
  /* assume all variables are written in full */
  Xt_uid *restrict uids = NULL, (*partDescUIDs)[numClients] = (Xt_uid(*)[numClients]) clientDeco.uids;
  Xt_idxlist(*partDesc)[numClients] = (Xt_idxlist(*)[numClients]) clientDeco.lists;
  int *partSizes = NULL;
  allocUIDLookup(numClients, &uids, &partSizes);
  for (int varID = 0; varID < nVars; ++varID)
    {
      Xt_idxlist outPartList = getOutDstList(vlistID, varID, conf);
      buildDecoPresetXmapsCommon(outPartList, numClients, partDescUIDs[varID], partDesc[varID], partSizes, uids, collComm, conf);
    }
  Free(uids);
}
#endif

/* build inventory of written variables for hypothetical full-stream */
static struct streamMapping *
streamMappingSpeculativeNew(int vlistID, struct partDescPreset clientDeco, const struct cdiPioConf *conf)
{
  (void) conf;
  int nVars = vlistNvars(vlistID);
  /* varMap[i] == index of header if variable i is written to,
   * numLvlsW[i] == number of levels of variable i or 0 if not written
   */
  size_t numWrittenRecords = 0;
  int *restrict numLvlsW = Malloc((size_t) nVars * sizeof(numLvlsW[0]));
  for (int varID = 0; varID < nVars; ++varID)
    {
      numWrittenRecords += (size_t) (numLvlsW[varID] = zaxisInqSize(vlistInqVarZaxis(vlistID, varID)));
    }
  /* set number of levels for each variable written to full */
  struct streamMapping *result = Malloc(sizeof(*result) + numWrittenRecords * sizeof(result->writtenRecords[0])
                                        + (size_t) nVars * 2 * sizeof(result->varMap[0]));
  struct recordWrite *restrict writtenRecords = result->writtenRecords;
  result->varMap = (void *) (writtenRecords + numWrittenRecords);
  result->numLvlsW = result->varMap + nVars;
  const int *conversion = clientDeco.conversion;
  {
    size_t j = (size_t) -1;
    /* initialized to shut up gcc, loop logic ensures initialization
     * to occur before first use */
    size_t recordNumElem = 0;
    int lastVarID = -1;
    for (int varID = 0; varID < nVars; ++varID)
      {
        size_t numLvl = (size_t) (result->numLvlsW[varID] = numLvlsW[varID]);
        if (varID != lastVarID)
          {
            int varShape[3];
            cdiPioQueryVarDims(varShape, vlistID, varID);
            recordNumElem = (size_t) varShape[0] * (size_t) varShape[1];
            lastVarID = varID;
          }
        result->varMap[varID] = 1;
        size_t elemSize = cdiPioElemSizeInference((size_t) varID, conversion);
        size_t recordDataSize = recordNumElem * elemSize;
        for (size_t lvl = 0; lvl < numLvl; ++lvl)
          writtenRecords[++j] = (struct recordWrite){ .varID = varID, .level = (int) lvl, .dataSize = recordDataSize };
      }
  }
  result->numVars = nVars;
  result->numWrittenRecords = (int) numWrittenRecords;
  Free(numLvlsW);
  result->layout = NULL;
  return result;
}

/* denote what will be aggregated at a single process */
struct passPlan
{
  unsigned recordAggStart, recordAggEnd;
  int varStart, varEnd;
};

struct passDict
{
  int varID;
  unsigned recordStart, recordEnd;
};

static size_t buildPassVarDict(size_t collSize, const struct passPlan *passes, const struct recordWrite *restrict writtenRecords,
                               size_t *dictSize, struct passDict **dict);

static size_t planPasses(const struct streamMapping *mapping, const struct cdiPioConf *conf, size_t collSize,
                         struct passPlan (**passes_)[collSize]);

struct idxlistAndSize
{
  Xt_idxlist list;
  int listSize;
};

static struct idxlistAndSize dstListFromRecordRange(int vlistID, int varID, int myVarStart, int myVarEnd, size_t myRecordStart,
                                                    size_t myRecordEnd, const struct recordWrite *restrict writtenRecords,
                                                    size_t recordStart, size_t recordEnd, const struct cdiPioConf *conf);

static void
buildDecoPresetXmapsGrib(int vlistID, struct partDescPreset clientDeco, MPI_Comm collComm, int collRank,
                         const struct cdiPioConf *conf)
{
  struct streamMapping *syntheticMapping = streamMappingSpeculativeNew(vlistID, clientDeco, conf);
  size_t collSize = (size_t) commInqSizeColl();
  struct passPlan(*passes)[collSize] = NULL;
  size_t numPasses = planPasses(syntheticMapping, conf, collSize, &passes);
  struct passDict *varsInPass = NULL;
  size_t varsInPassSize = 0;
  size_t numClients = (size_t) numClients_;
  /* assume all variables are written in full */
  Xt_uid *restrict uids = NULL, (*partDescUIDs)[numClients] = (Xt_uid(*)[numClients]) clientDeco.uids;
  Xt_idxlist(*partDesc)[numClients] = (Xt_idxlist(*)[numClients]) clientDeco.lists;
  int *restrict partSizes = NULL;
  allocUIDLookup(numClients, &uids, &partSizes);
  struct recordWrite *restrict writtenRecords = syntheticMapping->writtenRecords;
  for (size_t pass = 0; pass < numPasses; ++pass)
    {
      size_t myRecordStart = passes[pass][collRank].recordAggStart, myRecordEnd = passes[pass][collRank].recordAggEnd;
      size_t numVarsInPass = buildPassVarDict(collSize, passes[pass], writtenRecords, &varsInPassSize, &varsInPass);
      int myVarStart = passes[pass][collRank].varStart, myVarEnd = passes[pass][collRank].varEnd;
      for (size_t varIdx = 0; varIdx < numVarsInPass; ++varIdx)
        {
          int varID = varsInPass[varIdx].varID;
          struct idxlistAndSize dst
              = dstListFromRecordRange(vlistID, varID, myVarStart, myVarEnd, myRecordStart, myRecordEnd, writtenRecords,
                                       varsInPass[varIdx].recordStart, varsInPass[varIdx].recordEnd, conf);
          buildDecoPresetXmapsCommon(dst.list, numClients, partDescUIDs[varID], partDesc[varID], partSizes, uids, collComm, conf);
        }
    }
  Free(varsInPass);
  Free(uids);
  Free(syntheticMapping);
  Free(passes);
}

static void
buildDecoPresetXmaps(int streamID, struct partDescPreset clientDeco, MPI_Comm collComm, const struct cdiPioConf *conf)
{
  int filetype = streamInqFiletype(streamID);
  int vlistID = streamInqVlist(streamID);
  switch (filetype)
    {
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2:
      buildDecoPresetXmapsGrib(vlistID, clientDeco, collComm, commInqRankColl(), conf);
      /* writeGribStream(streamIdx, map, &data, &currentDataBufSize, conf); */
      break;
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
      {
        int rankOpen = cdiPioStream2Owner(streamID);
        outDstListConstruct createOutList;
#ifdef HAVE_PARALLEL_NC4
        if (rankOpen == CDI_PIO_COLLECTIVE_OPEN)
          createOutList = createOutDstListNetCDFPar;
        else
#endif
          {
            const int collRank = commInqRankColl();
            createOutList = rankOpen == collRank ? createOutDstListNetCDFSerialFunnel : createOutDstListNetCDFSerialAssist;
          }
        buildDecoPresetXmapsNetCDF(vlistID, clientDeco, collComm, createOutList, conf);
      }
#endif
    }
}

static Xt_redist
buildVarRedist(int headerIdx, size_t streamIdx,
               /* index list representing the data elements gathered on
                * this rank */
               Xt_idxlist dstList, Xt_idxlist *partDescPreset, const struct cdiPioConf *conf)
{
  const struct clientBuf *restrict clientBuf = rxWin[streamIdx].clientBuf;
  const struct winHeaderEntry *winDict = wHECast(clientBuf[0].mem);
  int varID = winDict[headerIdx].specific.dataRecord.varID;
  int conversion = winDict[headerIdx].id;
  size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
  size_t numClients = (size_t) numClients_;
  struct Xt_offset_ext *partExts = Malloc(numClients * sizeof(partExts[0]));
  MPI_Comm pioInterComm = cdiPioInqInterComm(), collComm = commInqCommColl();
  Xt_xmap gatherXmap;
  struct Xt_offset_ext gatherExt = { .start = 0, .size = xt_idxlist_get_num_indices(dstList), .stride = 1 };
  Xt_uid *restrict uids = NULL;
  int *restrict partSizes = NULL;
  bool cacheXmaps = conf->cacheXmaps;
  if (cacheXmaps)
    {
      allocUIDLookup(numClients, &uids, &partSizes);
      uids[0] = xt_idxlist_get_uid(dstList);
      for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx)
        {
          unsigned char *clientMem = clientBuf[clientIdx].mem;
          struct winHeaderEntry *partHeader = (wHECast(clientMem)) + headerIdx + 1;
          xassert(partHeader->id == PARTDESCMARKER);
          uids[clientIdx + 1] = unpackXTUID(partHeader->specific.partDesc.packedUID);
        }
      if ((gatherXmap = cdiPioXmapCacheLookup(XmapCache, uids, partSizes)))
        {
          for (size_t clientIdx = 0; clientIdx < numClients; ++clientIdx)
            {
              unsigned char *clientMem = clientBuf[clientIdx].mem;
              struct dataRecord *dataHeader = &(wHECast(clientMem))[headerIdx].specific.dataRecord;
              xassert(dataHeader->varID == varID && (wHECast(clientMem))[headerIdx].id == conversion);
              size_t charOfs = (size_t) ((clientMem + (wHECast(clientMem))[headerIdx].offset) - clientBuf[0].mem);
              int partSize = partSizes[clientIdx + 1];
              xassert(charOfs % elemSize == 0 && charOfs / elemSize + (size_t) partSize <= INT_MAX);
              int elemOfs = (int) (charOfs / elemSize);
              partExts[clientIdx].start = elemOfs;
              partExts[clientIdx].size = partSize;
              partExts[clientIdx].stride = 1;
            }
          goto finishXmapCaching;
        }
    }
  gatherXmap = buildVarXmap(partExts, clientBuf, (size_t) headerIdx, dstList, partDescPreset, pioInterComm, collComm, varID, conf);
  if (cacheXmaps)
    {
      partSizes[0] = gatherExt.size;
      for (size_t i = 0; i < numClients; ++i) partSizes[i + 1] = partExts[i].size;
      cdiPioXmapCacheAdd(XmapCache, uids, partSizes, gatherXmap);
    finishXmapCaching:
      Free(uids);
    }
  MPI_Datatype elemDt = conversion == DATA_HEADER_FLOAT ? MPI_FLOAT : MPI_DOUBLE;
  Xt_redist varRedist = xt_redist_p2p_ext_new(gatherXmap, (int) numClients, partExts, 1, &gatherExt, elemDt);
  if (!cacheXmaps) xt_xmap_delete(gatherXmap);
  Free(partExts);
  return varRedist;
}

static Xt_idxlist
buildVarSlicesIdxList(int vlistID, int varID, int startLvl, int numLvl, const struct cdiPioConf *conf)
{
  int varShape[3] = { 0, 0, 0 };
  cdiPioQueryVarDims(varShape, vlistID, varID);
  /* int varSize = varShape[0] * varShape[1] * varShape[2]; */
  Xt_int varShapeXt[3], origin[3] = { startLvl >= 0 ? (Xt_int) startLvl : 0, 0, 0 };
  int sliceShape[3];
  for (unsigned i = 0; i < 3; ++i) varShapeXt[2 - i] = (Xt_int) varShape[i];
  sliceShape[0] = numLvl >= 0 ? numLvl : (int) varShape[2];
  sliceShape[1] = varShape[1];
  sliceShape[2] = varShape[0];
  Xt_idxlist idxlist;
  if (conf->stripify)
    {
      /* FIXME: does not support grids larger than INT_MAX cells */
      int sliceSize = varShape[0] * varShape[1];
      Xt_int start = sliceSize * origin[0];
      int nstrides = sliceSize * sliceShape[0];
      idxlist = cdiPioIdxlistCacheAddStripes1(DstIdxlistCache, start, nstrides);
    }
  else
    idxlist = cdiPioIdxlistCacheAddSection3D(DstIdxlistCache, varShapeXt, origin, sliceShape);
  return idxlist;
}

static inline size_t
countMemMissingDouble(size_t n, const double *restrict data, double missVal)
{
  size_t nmiss = 0;
  for (size_t i = 0; i < n; ++i) nmiss += (data[i] == missVal);
  return nmiss;
}

static inline size_t
countMemMissingFloat(size_t n, const float *restrict data, double missVal)
{
  size_t nmiss = 0;
  float missValF = (float) missVal;
  if (missValF == missVal)
    for (size_t i = 0; i < n; ++i) nmiss += (data[i] == missValF);
  return nmiss;
}

static size_t
countVarChunkMissingVals(int vlistID, int varID, struct streamMapping *mapping, size_t chunkLen, int conversion,
                         const void *restrict data)
{
  size_t nmiss = 0;
  if (mapping->hasMissing[varID])
    {
      double missVal = vlistInqVarMissval(vlistID, varID);
      if (conversion == DATA_HEADER_DOUBLE)
        nmiss = countMemMissingDouble(chunkLen, data, missVal);
      else
        nmiss = countMemMissingFloat(chunkLen, data, missVal);
    }
  return nmiss;
}

static inline void
destructRetained(struct cacheRedist *restrict retained, size_t numRetained)
{
  for (size_t i = 0; i < (size_t) numRetained; ++i)
    if (retained[i].redist) xt_redist_delete(retained[i].redist);
}

/* return true if the new mapping for stream streamIdx is unchanged
 * with respect to the previous mapping such that a previously
 * constructed set of redists can be re-used.
 *
 * vlistID is only passed to prevent repeated calls to
 * streamInqVlist(streamID)
 */
static inline bool
handleRedistCache(size_t streamIdx, struct streamMapping *restrict mapping, size_t numPasses, int vlistID, MPI_Comm collComm)
{
  bool reuseRedists = false;
  if (!rxWin[streamIdx].retained)
    {
      rxWin[streamIdx].retained = Calloc(numPasses, sizeof(*rxWin[streamIdx].retained));
      rxWin[streamIdx].numRetained = numPasses;
      rxWin[streamIdx].prevLayout = mapping->layout;
      mapping->layout = NULL;
    }
  else
    {
      size_t numClients = (size_t) numClients_, numVars = (size_t) vlistNvars(vlistID);
      reuseRedists = !memcmp(mapping->layout, rxWin[streamIdx].prevLayout, numClients * numVars * sizeof(mapping->layout[0]));
      if (!reuseRedists)
        {
          Free(rxWin[streamIdx].prevLayout);
          rxWin[streamIdx].prevLayout = mapping->layout;
          mapping->layout = NULL;
        }
      {
        int temp = reuseRedists;
        xmpi(MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_INT, MPI_LAND, collComm));
        reuseRedists = temp;
      }
      if (!reuseRedists)
        {
          destructRetained(rxWin[streamIdx].retained, rxWin[streamIdx].numRetained);
          rxWin[streamIdx].retained = Realloc(rxWin[streamIdx].retained, numPasses * sizeof(*rxWin[streamIdx].retained));
          for (size_t i = 0; i < numPasses; ++i) rxWin[streamIdx].retained[i].redist = NULL;
          rxWin[streamIdx].numRetained = numPasses;
        }
    }
  return reuseRedists;
}

#ifdef HAVE_PARALLEL_NC4

#include <core/ppm_combinatorics.h>

/* collective writing variant */
static void
writeNetCDFStreamParallel(size_t streamIdx, struct streamMapping *mapping, void **data_, size_t *currentDataBufSize,
                          const struct cdiPioConf *conf)
{
  const int nvars = mapping->numVars;
  const int *restrict varMap = mapping->varMap;
  const int streamID = openStreams.entries[streamIdx], vlistID = streamInqVlist(streamID);
  const MPI_Comm collComm = commInqCommColl();
  // init redistCache if applicable
  const bool reuseRedists = conf->cacheRedists ? handleRedistCache(streamIdx, mapping, (size_t) nvars, vlistID, collComm) : false;
  struct cacheRedist *restrict retained = rxWin[streamIdx].retained;
  const struct clientBuf *restrict clientBuf = rxWin[streamIdx].clientBuf;
  const struct winHeaderEntry *winDict = wHECast(clientBuf[0].mem);
  Xt_idxlist *partDescPreset = rxWin[streamIdx].clientDeco.lists;

  for (int varID = 0; varID < nvars; ++varID)
    if (mapping->numLvlsW[varID])
      {
        int headerIdx = varMap[varID];
        size_t varSize;
        Xt_redist gatherRedist;
        struct PPM_extent varChunk[3];
        int myChunk[3][2];
        if (reuseRedists)
          {
            gatherRedist = retained[varID].redist;
            varSize = retained[varID].sliceSize;
            for (size_t i = 0; i < 3; ++i)
              {
                varChunk[i] = retained[varID].sliceExtent[i];
                myChunk[i][0] = PPM_extent_start(varChunk[i]);
                myChunk[i][1] = PPM_extent_end(varChunk[i]);
              }
          }
        else
          {
            /* prepare yaxt descriptor for write chunk */
            Xt_idxlist preWriteChunk;
            cdiPioNetCDFParChunk(vlistID, varID, &preWriteChunk, varChunk);
            for (int i = 0; i < 3; ++i)
              {
                myChunk[i][0] = PPM_extent_start(varChunk[i]);
                myChunk[i][1] = PPM_extent_end(varChunk[i]);
              }
            xdebug("Writing chunk { { %d, %d }, { %d, %d },"
                   " { %d, %d } }",
                   myChunk[0][0], myChunk[0][1], myChunk[1][0], myChunk[1][1], myChunk[2][0], myChunk[2][1]);
            varSize = (size_t) xt_idxlist_get_num_indices(preWriteChunk);
            /* transpose data into write deco */
            {
              Xt_idxlist *varPartDescPreset = partDescPreset ? partDescPreset + numClients_ * varID : NULL;
              gatherRedist = buildVarRedist(headerIdx, streamIdx, preWriteChunk, varPartDescPreset, conf);
              if (conf->cacheRedists)
                {
                  retained[varID].redist = gatherRedist;
                  retained[varID].sliceSize = varSize;
                  retained[varID].sliceExtent[0] = varChunk[0];
                  retained[varID].sliceExtent[1] = varChunk[1];
                  retained[varID].sliceExtent[2] = varChunk[2];
                }
            }
          }
        int conversion = winDict[headerIdx].id;
        size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
        resizeVarGatherBuf(elemSize * varSize, data_, currentDataBufSize);
        void *restrict data = *data_;
        xt_redist_s_exchange1(gatherRedist, rxWin[streamIdx].clientBuf[0].mem, data);

        if (!conf->cacheRedists) xt_redist_delete(gatherRedist);

        {
          /* count missing values if appropriate */
          size_t nmiss
              = countVarChunkMissingVals(vlistID, varID, mapping, (size_t) (PPM_extents_size(3, varChunk)), conversion, data);
          /* write chunk */
          if (conversion == DATA_HEADER_DOUBLE)
            streamWriteVarChunk(streamID, varID, (const int(*)[2]) myChunk, data, (int) nmiss);
          else
            streamWriteVarChunkF(streamID, varID, (const int(*)[2]) myChunk, data, (int) nmiss);
        }
      }
}

#endif

#if defined(HAVE_LIBNETCDF)
/* needed for writing when some files are only written to by a single process */
/* cdiOpenStreamMap(streamID) returns the writer process rank, or
 * CDI_PIO_COLLECTIVE if the file is written collectively */
int
cdiPioStream2Owner(int streamID)
{
  size_t streamIdx = indexOfID(&openStreams, streamID);
  xassert(streamIdx < SIZE_MAX);
  return rxWin[streamIdx].ownerRank;
}

/* for load-balancing purposes, count number of files per process */
/* cdiOpenFileCounts[rank] gives number of open files rank has to himself */
static int *cdiSerialOpenFileCount;

static int
cdiPioNextOpenRank(void)
{
  xassert(cdiSerialOpenFileCount != NULL);
  int commCollSize = commInqSizeColl();
  int minRank = 0, minOpenCount = cdiSerialOpenFileCount[0];
  for (int i = 1; i < commCollSize; ++i)
    if (cdiSerialOpenFileCount[i] < minOpenCount)
      {
        minOpenCount = cdiSerialOpenFileCount[i];
        minRank = i;
      }
  return minRank;
}

static void
cdiPioOpenFileOnRank(int rank)
{
  xassert(cdiSerialOpenFileCount != NULL && (unsigned) rank < (unsigned) commInqSizeColl());
  ++(cdiSerialOpenFileCount[rank]);
}

static void
cdiPioCloseFileOnRank(int rank)
{
  xassert(cdiSerialOpenFileCount != NULL && rank >= 0 && rank < commInqSizeColl());
  xassert(cdiSerialOpenFileCount[rank] > 0);
  --(cdiSerialOpenFileCount[rank]);
}

static void
cdiPioServerCdfDefVars(stream_t *streamptr)
{
  int rankOpen = cdiPioStream2Owner(streamptr->self);
  if (commInqIOMode() == PIO_NONE
#ifdef HAVE_PARALLEL_NC4
      || rankOpen == CDI_PIO_COLLECTIVE_OPEN
#endif
      || commInqRankColl() == rankOpen)
    cdfDefCoordinateVars(streamptr);
}

static void
writeNetCDFStreamSerial(size_t streamIdx, struct streamMapping *mapping, void **data_, size_t *currentDataBufSize,
                        const struct cdiPioConf *conf)
{
  const int nvars = mapping->numVars;
  const int *restrict varMap = mapping->varMap, *restrict numLvlsW = mapping->numLvlsW;
  /* determine process which has stream open (writer) and
   * which has data for which variable (var owner)
   * three cases need to be distinguished */
  const int streamID = openStreams.entries[streamIdx], vlistID = streamInqVlist(streamID);
  const int writerRank = cdiPioStream2Owner(streamID);
  const int collRank = commInqRankColl();
  const MPI_Comm collComm = commInqCommColl();  // HB: which communicator is to be supplied here?
  const bool reuseRedists
      = conf->cacheRedists != 0 ? handleRedistCache(streamIdx, mapping, (size_t) nvars, vlistID, collComm) : false;
  struct cacheRedist *restrict retained = rxWin[streamIdx].retained;

  const struct clientBuf *restrict clientBuf = rxWin[streamIdx].clientBuf;
  const struct winHeaderEntry *winDict = wHECast(clientBuf[0].mem);
  Xt_idxlist *partDescPreset = rxWin[streamIdx].clientDeco.lists;

  for (int varID = 0; varID < nvars; ++varID)
    if (numLvlsW[varID])
      {
        size_t varSize;
        Xt_redist gatherRedist;
        int headerIdx = varMap[varID];
        int conversion = winDict[headerIdx].id;
        size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);

        if (reuseRedists)
          {
            gatherRedist = retained[varID].redist;
            varSize = retained[varID].sliceSize;
            if (writerRank == collRank) resizeVarGatherBuf(varSize * elemSize, data_, currentDataBufSize);
          }
        else
          {
            Xt_idxlist dstList;
            if (writerRank == collRank)
              {
                dstList = buildVarSlicesIdxList(vlistID, varID, -1, -1, conf);
                varSize = (size_t) xt_idxlist_get_num_indices(dstList);
                resizeVarGatherBuf(varSize * elemSize, data_, currentDataBufSize);
              }
            else
              {
                varSize = 0;
                dstList = xt_idxempty_new();
              }
            Xt_idxlist *varPartDescPreset = partDescPreset ? partDescPreset + numClients_ * varID : NULL;
            gatherRedist = buildVarRedist(headerIdx, streamIdx, dstList, varPartDescPreset, conf);
            if (conf->cacheRedists)
              {
                retained[varID].redist = gatherRedist;
                retained[varID].sliceSize = varSize;
              }
          }
        void *restrict data = *data_;
        xt_redist_s_exchange1(gatherRedist, rxWin[streamIdx].clientBuf[0].mem, data);
        if (!conf->cacheRedists) xt_redist_delete(gatherRedist);
        if (writerRank == collRank)
          {
            size_t nmiss = countVarChunkMissingVals(vlistID, varID, mapping, varSize, conversion, data);
            if (conversion == DATA_HEADER_DOUBLE)
              streamWriteVar(streamID, varID, data, (int) nmiss);
            else
              streamWriteVarF(streamID, varID, data, (int) nmiss);
          }
      }
}

static void
writeNetCDFStream(size_t streamIdx, struct streamMapping *mapping, void **data_, size_t *currentDataBufSize,
                  const struct cdiPioConf *conf)
{
  void (*writeNetCDFStream_)(size_t streamIdx, struct streamMapping * mapping, void **data_, size_t *currentDataBufSize,
                             const struct cdiPioConf *conf)
      = writeNetCDFStreamSerial;
#ifdef HAVE_PARALLEL_NC4
  int streamID = openStreams.entries[streamIdx];
  int rankOpen = cdiPioStream2Owner(streamID);
  if (rankOpen == CDI_PIO_COLLECTIVE_OPEN) writeNetCDFStream_ = writeNetCDFStreamParallel;
#endif
  writeNetCDFStream_(streamIdx, mapping, data_, currentDataBufSize, conf);
}
#endif

static inline struct winHeaderEntry *
winDictEntry(size_t streamIdx, size_t client, size_t entry)
{
  return (wHECast(rxWin[streamIdx].clientBuf[client].mem)) + entry;
}

static struct streamMemLayout *
getLayout(size_t streamIdx)
{
  int streamID = openStreams.entries[streamIdx];
  size_t numClients = (size_t) numClients_;
  int vlistID = streamInqVlist(streamID);
  size_t numVars = (size_t) vlistNvars(vlistID);
  struct streamMemLayout(*layout)[numVars] = Calloc(numClients * numVars, sizeof(layout[0]));
  size_t numDataEntries = (size_t) (winDictEntry(streamIdx, 0, 0)->specific.headerSize.numDataEntries);
  for (size_t client = 0; client < numClients; ++client)
    for (size_t headerIdx = 1; headerIdx < numDataEntries; headerIdx += 2)
      {
        struct winHeaderEntry *varHeader = winDictEntry(streamIdx, client, headerIdx);
        xassert(varHeader->id == DATA_HEADER_DOUBLE || varHeader->id == DATA_HEADER_FLOAT);
        int conversion = varHeader->id;
        size_t varID = (size_t) varHeader[0].specific.dataRecord.varID;
        int offset = varHeader[0].offset;
        Xt_uid uid = unpackXTUID(varHeader[1].specific.partDesc.packedUID);
        layout[client][varID] = (struct streamMemLayout){ .varPartIdxListUID = uid, .offset = offset, .conversion = conversion };
      }
  return *layout;
}

/* build inventory of written variables for stream */
static struct streamMapping *
streamMappingNew(size_t streamIdx, const struct winHeaderEntry *winDict, const struct cdiPioConf *conf)
{
  int streamID = openStreams.entries[streamIdx];
  int numDataEntries = winDict[0].specific.headerSize.numDataEntries;
  int vlistID = streamInqVlist(streamID);
  int numVars = vlistNvars(vlistID);
  /* varMap[i] == index of header if variable i is written to,
   * numLvlsW[i] == number of levels of variable i or 0 if not written
   */
  int *restrict varMap = Calloc((size_t) numVars * 4, sizeof(varMap[0])), *restrict hasMissing = varMap + numVars,
                *restrict numLvlsW = varMap + 2 * numVars, *restrict hasMissing_ = varMap + 3 * numVars;
  for (int headerIdx = 1; headerIdx < numDataEntries; headerIdx += 2)
    {
      int varID = winDict[headerIdx].specific.dataRecord.varID;
      /* ensure a variable has not been enqueued twice */
      /* FIXME: this could better be ensured on client */
      xassert(varID < numVars && varID >= 0 && varMap[varID] == 0);
      varMap[varID] = headerIdx;
      hasMissing[varID] += winDict[headerIdx].specific.dataRecord.nmiss;
    }
  /* set numLvlsW[i] to 1 if varMap[i] != 0 on any collector,
   * also sets hasMissing_[i] to global reduction of hasMissing[i] */
  xmpi(MPI_Allreduce(varMap, numLvlsW, 2 * numVars, MPI_INT, MPI_LOR, commInqCommColl()));
  /* now find numbers of levels for each variable written anywhere */
  size_t numWrittenRecords = 0;
  for (int varID = 0; varID < numVars; ++varID)
    if (numLvlsW[varID]) numWrittenRecords += (size_t) (numLvlsW[varID] = zaxisInqSize(vlistInqVarZaxis(vlistID, varID)));
  struct streamMapping *result = Malloc(sizeof(*result) + numWrittenRecords * sizeof(result->writtenRecords[0])
                                        + (size_t) numVars * 3 * sizeof(result->varMap[0]));
  result->varMap = (void *) ((unsigned char *) result + sizeof(*result) + numWrittenRecords * sizeof(result->writtenRecords[0]));
  result->numLvlsW = result->varMap + numVars;
  result->hasMissing = result->varMap + 2 * numVars;
  {
    size_t j = (size_t) -1;
    /* initialized to shut up gcc, loop logic ensures initialization
     * to occur before first use */
    size_t recordNumElem = 0;
    int lastVarID = -1;
    for (int varID = 0; varID < numVars; ++varID)
      {
        size_t numLvl = (size_t) (result->numLvlsW[varID] = numLvlsW[varID]);
        if (varID != lastVarID)
          {
            int varShape[3];
            cdiPioQueryVarDims(varShape, vlistID, varID);
            recordNumElem = (size_t) varShape[0] * (size_t) varShape[1];
            lastVarID = varID;
          }
        size_t headerIdx = (size_t) varMap[varID];
        int conversion = winDict[headerIdx].id;
        size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
        size_t recordDataSize = recordNumElem * elemSize;
        result->varMap[varID] = varMap[varID];
        result->hasMissing[varID] = hasMissing_[varID];
        for (size_t lvl = 0; lvl < numLvl; ++lvl)
          result->writtenRecords[++j] = (struct recordWrite){ .varID = varID, .level = (int) lvl, .dataSize = recordDataSize };
      }
  }
  result->numVars = numVars;
  result->numWrittenRecords = (int) numWrittenRecords;
  Free(varMap);
  result->layout = conf->cacheRedists ? getLayout(streamIdx) : NULL;
  return result;
}

static void
streamMappingDelete(struct streamMapping **mapping)
{
  Free((*mapping)->layout);
  Free(*mapping);
  *mapping = NULL;
}

#if __GNUC__ == 11 && __GNUC_MINOR__ <= 2
/* gcc 11.1 has a bug in the -fsanitize=undefined functionality
 * which creates a bogus warning without the below suppression, bug
 * report at https://gcc.gnu.org/bugzilla/show_bug.cgi?id=101585 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla-parameter"
#endif

/**
 * @param[out] passes_ pointer to pointer to 2-dimensional array of
 * records of dimensions $number of passes \cdot number of collectors$,
 * where $(*passes_)[pass][i]$ details the records written by collector
 * rank \a i
 * @return number of passes
 */
static size_t
planPasses(const struct streamMapping *mapping, const struct cdiPioConf *conf, size_t collSize,
           struct passPlan (**passes_)[collSize])
{
  size_t numPasses = 0;
  size_t recordAggBufLim = conf->recordAggBufLimMB * 1024 * 1024, totalAggBufSpace = recordAggBufLim * collSize, totalWritten = 0;
  /* find total size of data written for the stream and build prefix sums */
  size_t numWrittenRecords = (size_t) mapping->numWrittenRecords;

  if (numWrittenRecords == 0) return 0;
  size_t *restrict recordDataSizePfxSums = Malloc((numWrittenRecords + 1 + collSize + 1) * sizeof(*recordDataSizePfxSums)),
                   *restrict recordSeparations = recordDataSizePfxSums + numWrittenRecords + 1;
  const struct recordWrite *restrict writtenRecords = mapping->writtenRecords;

  recordDataSizePfxSums[0] = 0;
  for (size_t i = 0; i < numWrittenRecords; ++i)
    {
      size_t recordDataSize = writtenRecords[i].dataSize;
      recordDataSizePfxSums[i + 1] = recordDataSizePfxSums[i] + recordDataSize;
      totalWritten += recordDataSize;
    }
  /* move if into loop for handling last pass */
  if (totalWritten < totalAggBufSpace)
    {
      /* don't go to limit of some tasks where a single pass will be
       * sufficient to write everything, compute load-balancing
       * instead */
      numPasses = 1;
      struct passPlan *passes = Malloc(sizeof(*passes) * collSize);
      cdiPioDeco1D_CCP(numWrittenRecords, recordDataSizePfxSums, collSize, recordSeparations);
      for (size_t rank = 0; rank < collSize; ++rank)
        {
          size_t startRecord = recordSeparations[rank], lastRecord = recordSeparations[rank + 1] - 1;
          passes[rank] = (struct passPlan){
            .recordAggStart = (unsigned) startRecord,
            .recordAggEnd = (unsigned) lastRecord,
            .varStart = writtenRecords[startRecord].varID,
            .varEnd = writtenRecords[lastRecord].varID,
          };
        }
      *passes_ = (struct passPlan(*)[collSize]) passes;
    }
  else
    {
      /* aggregate as many records on each task to fill up to
       * recordAggLim data bytes, but use at least one, unless none
       * remain */
      size_t firstRecordOfPass = 0, curRecord;
#if !defined __PGI || __PGIC__ != 20
      struct passPlan(*passes)[collSize] = NULL;
#else
      /* we need the following workaround because pggpp2 fails to produce
       * LLVM IR assembly file with -O2 optimization level */
      struct passPlan *passes = NULL;
#endif
      size_t sizeof_pass = collSize * sizeof(struct passPlan);
      do
        {
          size_t taskBegin = firstRecordOfPass;
          curRecord = firstRecordOfPass - 1;
          passes = Realloc(passes, sizeof_pass * (numPasses + 1));
          for (size_t rank = 0; rank < collSize; ++rank)
            {
              size_t recordAggBufSize = 0;
              while (curRecord + 1 < numWrittenRecords
                     && ((recordAggBufSize + writtenRecords[curRecord + 1].dataSize) < recordAggBufLim))
                recordAggBufSize += writtenRecords[++curRecord].dataSize;
              if (curRecord == taskBegin - 1 && curRecord + 1 < numWrittenRecords) ++curRecord;
#if !defined __PGI || __PGIC__ != 20
              passes[numPasses][rank]
#else
              passes[numPasses * collSize + rank]
#endif
                  = (struct passPlan){
                      .recordAggStart = (unsigned) taskBegin,
                      .recordAggEnd = (unsigned) curRecord,
                      .varStart = writtenRecords[taskBegin].varID,
                      .varEnd = writtenRecords[curRecord].varID,
                    };
              taskBegin = curRecord + 1;
            }
          ++numPasses, firstRecordOfPass = curRecord + 1;
        }
      while (curRecord + 1 < numWrittenRecords);
#if !defined __PGI || __PGIC__ != 20
      *passes_ = passes;
#else
      *passes_ = (struct passPlan(*)[collSize]) passes;
#endif
    }
  Free(recordDataSizePfxSums);
  return numPasses;
}

#if __GNUC__ == 11 && __GNUC_MINOR__ <= 2
#pragma GCC diagnostic pop
#endif

static inline size_t
szmin(size_t a, size_t b)
{
  return a <= b ? a : b;
}

static inline size_t
szmax(size_t a, size_t b)
{
  return a >= b ? a : b;
}

static size_t
aggBufAppend(int fileID, const void *restrict ptr, size_t size)
{
  (void) fileID;
  size_t aggBufSize = aggBuf.size, aggBufUsed = aggBuf.used;
  void *restrict buf = aggBuf.mem;
  if (aggBufUsed + size > aggBufSize) aggBuf.mem = buf = Realloc(buf, (aggBuf.size = aggBufUsed + size));
  memcpy((unsigned char *) buf + aggBufUsed, ptr, size);
  aggBuf.used = aggBufUsed + size;
  return size;
}

static void
aggBufFlush(int streamID, int fileID, size_t (*cdiPioFileWrite)(int, const void *restrict, size_t, int))
{
  cdiPioFileWrite(fileID, aggBuf.mem, aggBuf.used, streamInqCurTimestepID(streamID));
  aggBuf.used = 0;
}

static size_t
buildPassVarDict(size_t collSize, const struct passPlan *passes, const struct recordWrite *restrict writtenRecords,
                 size_t *dictSize, struct passDict **dict)
{
  unsigned base = passes[0].recordAggStart;
  size_t numRecordsInPass = passes[collSize - 1].recordAggEnd - base + 1;
  size_t maxVarsInPass = (size_t) (passes[collSize - 1].varEnd - passes[0].varStart + 1);
  struct passDict *varsInPass = *dict;
  size_t oldSize = *dictSize, newSize = szmin(numRecordsInPass, maxVarsInPass);
  if (oldSize < newSize)
    {
      *dictSize = newSize;
      varsInPass = *dict = Realloc(varsInPass, newSize * sizeof(*varsInPass));
    }
  /* establish variables involved in this pass */
  size_t numVarsInPass = 1;
  varsInPass[0].recordStart = base;
  int lastSeenVarID = varsInPass[0].varID = writtenRecords[base].varID;
  for (size_t i = 1; i < numRecordsInPass; ++i)
    if (lastSeenVarID != writtenRecords[base + i].varID)
      {
        varsInPass[numVarsInPass - 1].recordEnd = (unsigned) (base + i - 1);
        varsInPass[numVarsInPass].varID = lastSeenVarID = writtenRecords[base + i].varID;
        varsInPass[numVarsInPass].recordStart = (unsigned) (base + i);
        ++numVarsInPass;
      }
  varsInPass[numVarsInPass - 1].recordEnd = (unsigned) (base + numRecordsInPass - 1);
  return numVarsInPass;
}

static struct idxlistAndSize
dstListFromRecordRange(int vlistID, int varID, int myVarStart, int myVarEnd, size_t myRecordStart, size_t myRecordEnd,
                       const struct recordWrite *restrict writtenRecords, size_t recordStart, size_t recordEnd,
                       const struct cdiPioConf *conf)
{
  /* is this process writing part of this variable? */
  Xt_idxlist dstList;
  int listSize;
  if (myRecordStart <= myRecordEnd && myVarStart <= varID && myVarEnd >= varID)
    {
      size_t myVarRecordStart = writtenRecords[myRecordStart].varID == varID ? myRecordStart : recordStart;
      size_t myLevelStart = (size_t) writtenRecords[myVarRecordStart].level;
      size_t myVarRecordEnd = writtenRecords[myRecordEnd].varID == varID ? myRecordEnd : recordEnd;
      size_t myNumLevels = (size_t) writtenRecords[myVarRecordEnd].level - myLevelStart + 1;
      dstList = buildVarSlicesIdxList(vlistID, varID, (int) myLevelStart, (int) myNumLevels, conf);
      listSize = xt_idxlist_get_num_indices(dstList);
    }
  else
    {
      dstList = xt_idxempty_new();
      listSize = 0;
    }
  return (struct idxlistAndSize){ .list = dstList, .listSize = listSize };
}

static void
writeGribStream(size_t streamIdx, struct streamMapping *mapping, void **data_, size_t *currentDataBufSize,
                const struct cdiPioConf *conf)
{
  const struct clientBuf *restrict clientBuf = rxWin[streamIdx].clientBuf;
  int streamID = openStreams.entries[streamIdx];
  int vlistID = streamInqVlist(streamID);
  int fileID = streamInqFileID(streamID);
  MPI_Comm collComm = commInqCommColl();
  size_t collSize = (size_t) commInqSizeColl();
  size_t collRank = (size_t) commInqRankColl();
  struct passPlan(*passes)[collSize] = NULL;
  size_t numPasses = planPasses(mapping, conf, collSize, &passes);
  Xt_redist *varRedists = NULL;
  size_t numClients = (size_t) numClients_;
  Xt_idxlist(*partDescPreset)[numClients] = (Xt_idxlist(*)[numClients]) rxWin[streamIdx].clientDeco.lists;
  struct recordWrite *restrict writtenRecords = mapping->writtenRecords;
  size_t (*cdiPioFileWrite)(int fileID, const void *restrict buffer, size_t len, int tsID)
      = (size_t(*)(int, const void *restrict, size_t, int)) namespaceSwitchGet(NSSWITCH_FILE_WRITE).func;
  bool reuseRedists
      = conf->cacheRedists != 0 ? handleRedistCache(streamIdx, mapping, (size_t) numPasses, vlistID, collComm) : false;
  struct cacheRedist *restrict retained = rxWin[streamIdx].retained;
  struct passDict *varsInPass = NULL;
  size_t varsInPassSize = 0, maxNumVarsAlloc = 0;
  MPI_Aint *displ = NULL;
  const struct winHeaderEntry *winDict = wHECast(clientBuf[0].mem);
  for (size_t pass = 0; pass < numPasses; ++pass)
    {
      size_t myRecordStart = passes[pass][collRank].recordAggStart, myRecordEnd = passes[pass][collRank].recordAggEnd;
      size_t myAggSize = 0;
      /* build or fetch from cache redists for all variables involved in current write pass */
      Xt_redist compositePassRedist;
      if (reuseRedists)
        {
          compositePassRedist = retained[pass].redist;
          myAggSize = (size_t) retained[pass].sliceSize;
        }
      else
        {
          size_t numVarsInPass = buildPassVarDict(collSize, passes[pass], writtenRecords, &varsInPassSize, &varsInPass);
          int myVarStart = passes[pass][collRank].varStart, myVarEnd = passes[pass][collRank].varEnd;
          if (numVarsInPass > maxNumVarsAlloc)
            {
              maxNumVarsAlloc = numVarsInPass;
              varRedists = Realloc(varRedists, numVarsInPass * sizeof(*varRedists));
              displ = Realloc(displ, (numVarsInPass * 2 + 1) * sizeof(*displ));
            }
          memset(displ, 0, sizeof(*displ) * (numVarsInPass + 1));
          for (size_t varIdx = 0; varIdx < numVarsInPass; ++varIdx)
            {
              int varID = varsInPass[varIdx].varID, headerIdx = mapping->varMap[varID];
              int conversion = winDict[headerIdx].id;
              size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
              struct idxlistAndSize dst
                  = dstListFromRecordRange(vlistID, varID, myVarStart, myVarEnd, myRecordStart, myRecordEnd, writtenRecords,
                                           varsInPass[varIdx].recordStart, varsInPass[varIdx].recordEnd, conf);
              myAggSize += (size_t) dst.listSize * elemSize;
              displ[numVarsInPass + varIdx + 1] = (MPI_Aint) myAggSize;
              Xt_idxlist *varPartDescPreset = partDescPreset ? partDescPreset[varIdx] : NULL;
              varRedists[varIdx] = buildVarRedist(headerIdx, streamIdx, dst.list, varPartDescPreset, conf);
            }
          /* merge all redists for current pass */
          if (numVarsInPass > 1)
            {
              compositePassRedist
                  = xt_redist_collection_static_new(varRedists, (int) numVarsInPass, displ, displ + numVarsInPass, collComm);
              /* free individual redists */
              for (size_t varIdx = 0; varIdx < numVarsInPass; ++varIdx) xt_redist_delete(varRedists[varIdx]);
            }
          else
            compositePassRedist = varRedists[0];
          if (conf->cacheRedists)
            {
              retained[pass].redist = compositePassRedist;
              retained[pass].sliceSize = myAggSize;
            }
        }
      /* resize gather buffer if needed */
      resizeVarGatherBuf(myAggSize, data_, currentDataBufSize);
      /* execute composite redist */
      xt_redist_s_exchange1(compositePassRedist, clientBuf[0].mem, *data_);
      /* delete composite redist */
      if (!conf->cacheRedists) xt_redist_delete(compositePassRedist);

      /* append encoded data records from this pass to buffer written later */
      /* todo: develop better heuristic for buffer size */
      if (myAggSize > aggBuf.size)
        {
          Free(aggBuf.mem);
          size_t aggBufSize = szmax((size_t) conf->recordAggBufLimMB * (size_t) 1024 * (size_t) 1024, myAggSize);
          if (posix_memalign(&aggBuf.mem, cdiGetPageSize(conf->largePageAlign), aggBufSize) == 0)
            ;
          else
            aggBuf.mem = Malloc(aggBufSize);
          aggBuf.size = aggBufSize;
        }
      namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(aggBufAppend));
      /* write records to aggregation buffer */
      if (myRecordStart <= myRecordEnd)
        {
          size_t varIdx = (size_t) -1;
          int varID = -1;
          size_t recordDataOfs = 0;
          const unsigned char *restrict data = *data_;
          for (size_t recordIdx = myRecordStart; recordIdx <= myRecordEnd; ++recordIdx)
            {
              int level = writtenRecords[recordIdx].level;
              int prevVarID = varID;
              varID = writtenRecords[recordIdx].varID;
              varIdx += varID != prevVarID;
              size_t recordSize = writtenRecords[recordIdx].dataSize;
              int headerIdx = mapping->varMap[varID];
              int conversion = winDict[headerIdx].id;
              size_t elemSize = conversion == DATA_HEADER_FLOAT ? sizeof(float) : sizeof(double);
              size_t nvals = recordSize / elemSize;
              size_t nmiss = countVarChunkMissingVals(vlistID, varID, mapping, nvals, conversion, data + recordDataOfs);
              if (conversion == DATA_HEADER_DOUBLE)
                streamWriteVarSlice(streamID, varID, level, (const double *) (const void *) (data + recordDataOfs), (int) nmiss);
              else
                streamWriteVarSliceF(streamID, varID, level, (const float *) (const void *) (data + recordDataOfs), (int) nmiss);
              recordDataOfs += recordSize;
            }
          aggBufFlush(streamID, fileID, cdiPioFileWrite);
        }
      else
        /* write zero bytes to trigger synchronization code in fileWrite */
        cdiPioFileWrite(fileID, NULL, 0, streamInqCurTimestepID(streamID));
      namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(cdiPioFileWrite));
    }
  Free(displ);
  Free(varRedists);
  Free(varsInPass);
  Free(passes);
}

static void
readGetBuffers(size_t streamIdx, const struct cdiPioConf *conf)
{
  int streamID = openStreams.entries[streamIdx];
  xdebug("%s", "START");

  struct winHeaderEntry *winDict = wHECast(rxWin[streamIdx].clientBuf[0].mem);
  xassert(winDict[0].id == HEADERSIZEMARKER);
  {
    int dictSize = rxWin[streamIdx].clientBuf[0].dictSize,
        firstNonRPCEntry = dictSize - winDict[0].specific.headerSize.numRPCEntries - 1, headerIdx, numFuncCalls = 0;
    for (headerIdx = dictSize - 1; headerIdx > firstNonRPCEntry; --headerIdx)
      {
        xassert(winDict[headerIdx].id >= MINFUNCID && winDict[headerIdx].id <= MAXFUNCID);
        ++numFuncCalls;
        readFuncCall(winDict + headerIdx, streamIdx);
      }
    xassert(numFuncCalls == winDict[0].specific.headerSize.numRPCEntries);
  }
  /* build list of streams, data was transferred for */
  {
    struct streamMapping *map = streamMappingNew(streamIdx, winDict, conf);
    /* TODO: build list of rma buffer layout here to check if caching can be done */
    void *data = NULL;
    size_t currentDataBufSize = 0;
    int filetype = streamInqFiletype(streamID);

    switch (cdiBaseFiletype(filetype))
      {
      case CDI_FILETYPE_GRIB: writeGribStream(streamIdx, map, &data, &currentDataBufSize, conf); break;
#ifdef HAVE_LIBNETCDF
      case CDI_FILETYPE_NETCDF: writeNetCDFStream(streamIdx, map, &data, &currentDataBufSize, conf); break;
#endif
      default: xabort("unhandled filetype in parallel I/O.");
      }
    streamMappingDelete(&map);
    Free(map);
    Free(data);
  }
  xdebug("%s", "RETURN");
}

/************************************************************************/

static void
clearModelWinBuffer(size_t streamIdx)
{
  xassert(streamIdx < openStreams.size && rxWin != NULL && rxWin[streamIdx].clientBuf[0].mem != NULL);
  size_t clearSize;
#if 0
  /* for debugging purposes it might be smart to set the whole buffer
   * to 0 */
  clearSize = (size_t)(rxWin[streamIdx].clientBuf[numClients_ - 1].mem
                          - rxWin[streamIdx].clientBuf[0].mem)
    + rxWin[streamIdx].clientBuf[numClients_ - 1].size;
#else
  /* but normally, setting the first header record to zero will suffice */
  clearSize = sizeof(struct winHeaderEntry);
#endif
  memset(rxWin[streamIdx].clientBuf[0].mem, 0, clearSize);
}

/************************************************************************/

static void
getTimeStepData(int *streamActivity, const struct cdiPioConf *conf)
{
  MPI_Group clientGroup = cdiPioInqRemoteGroup();

  xdebug("%s", "START");

  for (size_t streamIdx = 0; streamIdx < openStreams.size; ++streamIdx)
    if (streamActivity[streamIdx])
      {
        clearModelWinBuffer(streamIdx);
        xmpi(MPI_Win_start(clientGroup, 0, rxWin[streamIdx].getWin));
        /* FIXME: this needs to use MPI_PACKED for portability */
        for (size_t i = 0; i < (size_t) numClients_; ++i)
          xmpi(MPI_Get(rxWin[streamIdx].clientBuf[i].mem, (int) rxWin[streamIdx].clientBuf[i].size, MPI_UNSIGNED_CHAR,
                       clientRanks_[i], 0, (int) rxWin[streamIdx].clientBuf[i].size, MPI_UNSIGNED_CHAR, rxWin[streamIdx].getWin));
        xmpi(MPI_Win_complete(rxWin[streamIdx].getWin));
        if (!conf->batchedRMA) readGetBuffers(streamIdx, conf);
      }
  if (conf->batchedRMA)
    for (size_t streamIdx = 0; streamIdx < openStreams.size; ++streamIdx)
      if (streamActivity[streamIdx]) readGetBuffers(streamIdx, conf);

  xdebug("%s", "RETURN");
}

/************************************************************************/

static int
cdiPioServerStreamOpen(const char *filename, char filemode, int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
{
  int fileID = -1;
#ifdef HAVE_LIBNETCDF
  /* Only needs initialization to shut up gcc */
  int rank = -1;
#endif
  switch (filetype)
    {
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
#ifdef HAVE_PARALLEL_NC4
      {
        struct cdiPioNcCreateLongJmpRetBuf retJmpBuf;
        retJmpBuf.openRank = cdiPioNextOpenRank();
#if !defined TLS && defined HAVE_PTHREAD
        pthread_setspecific(cdiPioCdfJmpKey, &jmpBuf);
#else
        cdiPioCdfJmpBuf = &retJmpBuf;
#endif
        if (!setjmp(retJmpBuf.jmpBuf)) /* attempt parallel open first */
          /* in case it fails, ranks other than retJmpBuf.openRank
           * will call longjmp and return 1 from the above setjmp */
          fileID = cdiStreamOpenDefaultDelegate(filename, filemode, filetype, streamptr, recordBufIsToBeCreated);
        rank = retJmpBuf.openRank;
        if (rank != CDI_PIO_COLLECTIVE_OPEN)
          {
            streamptr->filetype = filetype;
            if (commInqIOMode() != PIO_NONE) xmpi(MPI_Bcast(&fileID, 1, MPI_INT, rank, commInqCommColl()));
            cdiPioOpenFileOnRank(rank);
          }
      }
#else
      {
        int ioMode = commInqIOMode();
        if (ioMode == PIO_NONE || commInqRankColl() == (rank = cdiPioNextOpenRank()))
          fileID = cdiStreamOpenDefaultDelegate(filename, filemode, filetype, streamptr, recordBufIsToBeCreated);
        else
          streamptr->filetype = filetype;
        if (ioMode != PIO_NONE) xmpi(MPI_Bcast(&fileID, 1, MPI_INT, rank, commInqCommColl()));
        cdiPioOpenFileOnRank(rank);
      }
#endif
      break;
#endif
    default: fileID = cdiStreamOpenDefaultDelegate(filename, filemode, filetype, streamptr, recordBufIsToBeCreated);
    }
  if (fileID >= 0)
    {
      size_t oldNumStreams = openStreams.size;
      size_t streamIdx = insertID(&openStreams, streamptr->self);
      size_t numStreams = openStreams.size;
      struct clientBuf *oldClientBufs = rxWin ? rxWin[0].clientBuf : NULL;
      rxWin = Realloc(rxWin, numStreams * sizeof(rxWin[0]));
      struct clientBuf *restrict newClientBufs
          = Realloc(oldClientBufs, sizeof(rxWin[0].clientBuf[0]) * (size_t) numClients_ * numStreams);
      if (newClientBufs != oldClientBufs)
        for (size_t i = 0; i < numStreams; ++i) rxWin[i].clientBuf = newClientBufs + i * (size_t) numClients_;
      else if (oldNumStreams < numStreams)
        for (size_t i = oldNumStreams; i < numStreams; ++i) rxWin[i].clientBuf = newClientBufs + i * (size_t) numClients_;
      rxWin[streamIdx].getWin = MPI_WIN_NULL;
      rxWin[streamIdx].clientBuf[0].mem = NULL;
      rxWin[streamIdx].prevLayout = NULL;
      rxWin[streamIdx].retained = NULL;
      rxWin[streamIdx].numRetained = 0;
      rxWin[streamIdx].clientDeco.lists = NULL;
      rxWin[streamIdx].clientDeco.uids = NULL;
      rxWin[streamIdx].clientDeco.conversion = NULL;
#ifdef HAVE_LIBNETCDF
      rxWin[streamIdx].ownerRank = rank;
#endif
    }
  return fileID;
}

static size_t
getMaxNumStreamWrites(stream_t *streamptr)
{
  int filetype = streamptr->filetype;
  int vlistID = streamptr->vlistID;
  size_t numVars = (size_t) vlistNvars(streamptr->vlistID), maxNumStreamWrites = 0, numColl = (size_t) commInqSizeColl();
  switch (filetype)
    {
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
      {
        int rankOpen = cdiPioStream2Owner(streamptr->self);
        if (commInqIOMode() == PIO_NONE
#ifdef HAVE_PARALLEL_NC4
            || rankOpen == CDI_PIO_COLLECTIVE_OPEN
#endif
            || (commInqRankColl() == rankOpen))
          maxNumStreamWrites = (numVars + numColl - 1) / numColl;
      }
      break;
#endif
    default:
      for (size_t varID = 0; varID < numVars; ++varID)
        {
          size_t numRec = (size_t) (zaxisInqSize(vlistInqVarZaxis(vlistID, (int) varID)));
          maxNumStreamWrites += (numRec + numColl - 1) / numColl;
        }
    }
  return maxNumStreamWrites;
}

static void
cdiPioServerStreamClose(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  int fileID = streamptr->fileID;
  int filetype = streamptr->filetype;
  int vlistID = streamptr->vlistID;
  const struct cdiPioConf *conf = cdiPioGetConf();
  if (vlistID != CDI_UNDEFID)
    {
      size_t maxNumStreamWrites = getMaxNumStreamWrites(streamptr);
      neededDstIdxlistCacheSize -= maxNumStreamWrites;
      if (conf->cacheXmaps) neededXmapCacheSize -= maxNumStreamWrites;
    }
  if (fileID == CDI_UNDEFID)
    Warning("File %s not open!", streamptr->filename);
  else
    {
      switch (cdiBaseFiletype(filetype))
        {
#ifdef HAVE_LIBNETCDF
        case CDI_FILETYPE_NETCDF:
          {
            int rankOpen = cdiPioStream2Owner(streamptr->self);
            if (commInqIOMode() == PIO_NONE
#ifdef HAVE_PARALLEL_NC4
                || rankOpen == CDI_PIO_COLLECTIVE_OPEN
#endif
                || commInqRankColl() == rankOpen)
              cdiStreamCloseDefaultDelegate(streamptr, recordBufIsToBeDeleted);
#ifdef HAVE_PARALLEL_NC4
            if (rankOpen != CDI_PIO_COLLECTIVE_OPEN)
#endif
              cdiPioCloseFileOnRank(rankOpen);
          }
          break;
#endif
        default: cdiStreamCloseDefaultDelegate(streamptr, recordBufIsToBeDeleted);
        }
      int streamID = streamptr->self;
      size_t streamIdx = indexOfID(&openStreams, streamID);
      destructRetained(rxWin[streamIdx].retained, rxWin[streamIdx].numRetained);
      if (rxWin[streamIdx].clientDeco.lists)
        cdiPioDestroyPartDescPreset((size_t) numClients_, (size_t) (vlistNvars(vlistID)), &rxWin[streamIdx].clientDeco);
      Free(rxWin[streamIdx].retained);
      Free(rxWin[streamIdx].prevLayout);
      cdiPioServerStreamWinDestroy(streamIdx, conf);
      removeID(&openStreams, streamID);
    }
  void (*streamCloseCallBack)(int streamID) = (void (*)(int)) conf->callbacks[CDIPIO_CALLBACK_POSTSTREAMCLOSE];
  streamCloseCallBack(streamptr->self);
}

#if defined HAVE_H5GET_LIBVERSION && defined HAVE_PARALLEL_NC4
static bool parH5ZeroCountProblem = false;
#endif

#ifdef HAVE_LIBNETCDF
static void
cdiPioCdfDefTimestep(stream_t *streamptr, int tsID, size_t valCount)
{
  int streamID = streamptr->self, rankOpen = cdiPioStream2Owner(streamID);
#ifdef HAVE_PARALLEL_NC4
  valCount = rankOpen != CDI_PIO_COLLECTIVE_OPEN
#if defined HAVE_H5GET_LIBVERSION
             || (parH5ZeroCountProblem && streamptr->filetype == CDI_FILETYPE_NC4)
#else
             || streamptr->filetype == CDI_FILETYPE_NC4
#endif
             || commInqRankColl() == 0;
#endif
  if (commInqIOMode() == PIO_NONE
#ifdef HAVE_PARALLEL_NC4
      || rankOpen == CDI_PIO_COLLECTIVE_OPEN
#endif
      || commInqRankColl() == rankOpen)
    cdfDefTimestep(streamptr, tsID, valCount);
}
#endif

static void
cdiPioRecvStreamOpen(void *buffer, int size, int *pos, MPI_Comm pioInterComm)
{
  int clientStreamID, filetype, fname_len;
  {
    int soHdr[3];
    xmpi(MPI_Unpack(buffer, size, pos, soHdr, 3, MPI_INT, pioInterComm));
    clientStreamID = soHdr[0];
    filetype = soHdr[1];
    fname_len = soHdr[2];
  }
  char filemode;
  xmpi(MPI_Unpack(buffer, size, pos, &filemode, 1, MPI_CHAR, pioInterComm));
  MPI_Request *requests = Malloc((size_t) numClients_ * sizeof(requests[0]) + (size_t) fname_len + 1);
  char *filename = (char *) ((unsigned char *) requests + (size_t) numClients_ * sizeof(requests[0]));
  xmpi(MPI_Unpack(buffer, size, pos, filename, fname_len, MPI_CHAR, pioInterComm));
  filename[fname_len] = '\0';
  xassert(filemode == 'w');
  int serverStreamID = namespaceAdaptKey2(clientStreamID);
  int curStatus = reshGetStatus(serverStreamID, &streamOps);
  xassert(!(curStatus & RESH_IN_USE_BIT));
  int streamID = streamOpenID(filename, filemode, filetype, serverStreamID);
  int fileID = (streamID >= 0) ? streamInqFileID(streamID) : streamID;
  for (size_t i = 0; i < (size_t) numClients_; ++i)
    xmpi(MPI_Isend(&fileID, 1, MPI_INT, clientRanks_[i], STREAMOPEN, pioInterComm, requests + i));
  xmpi(MPI_Waitall(numClients_, requests, MPI_STATUSES_IGNORE));
  Free(requests);
}

static void
cdiPioRecvStreamClose(void *buffer, int size, int *pos, MPI_Comm pioInterComm, bool flushNeeded, const struct cdiPioConf *conf)
{
  int clientStreamID;
  xmpi(MPI_Unpack(buffer, size, pos, &clientStreamID, 1, MPI_INT, pioInterComm));
  int serverStreamID = namespaceAdaptKey2(clientStreamID);
  if (flushNeeded)
    {
      int streamActivity[openStreams.size];
      size_t streamIdx = indexOfID(&openStreams, serverStreamID);
      for (size_t i = 0; i < openStreams.size; ++i) streamActivity[i] = 0;
      streamActivity[streamIdx] = 1;
      getTimeStepData(streamActivity, conf);
      conf->callbacks[CDIPIO_CALLBACK_POSTWRITEBATCH]();
    }
  streamClose(serverStreamID);
}

#ifdef HAVE_LIBNETCDF
static void
cdiPioCdfGridAccess(int streamID, int vlistID)
{
  int filetype = streamInqFiletype(streamID);
  switch (filetype)
    {
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
#ifdef HAVE_PARALLEL_NC4
      if (cdiPioStream2Owner(streamID) != CDI_PIO_COLLECTIVE_OPEN)
#endif
        {
          int nGrids = vlistNgrids(vlistID);
          for (int gridIdx = 0; gridIdx < nGrids; ++gridIdx)
            {
              int gridID = vlistGrid(vlistID, gridIdx);
              if (!cdiPioDistGridIndividualQueriesEnabled(gridID)) cdiPioDistGridEnableIndividualQueries(gridID);
            }
        }
    }
}
#else
#define cdiPioCdfGridAccess(streamID, vlistID)
#endif

static void
cdiPioRecvStreamDefVlist(void *buffer, int size, int *pos, MPI_Comm pioInterComm, int tag, const struct cdiPioConf *conf)
{
  int serverStreamID, serverVlistID;
  {
    int msgData[defVlistNInts];
    xmpi(MPI_Unpack(buffer, size, pos, &msgData, defVlistNInts, MPI_INT, pioInterComm));
    serverStreamID = namespaceAdaptKey2(msgData[0]);
    serverVlistID = namespaceAdaptKey2(msgData[1]);
  }
  stream_t *streamptr = stream_to_pointer(serverStreamID);
  cdiPioCdfGridAccess(serverStreamID, serverVlistID);
  cdiStreamSetupVlist(streamptr, serverVlistID);
  size_t streamIdx = indexOfID(&openStreams, serverStreamID);
  cdiPioServerStreamWinCreate(streamIdx, cdiPioInqCollClientIntraComm());
  int numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl();
  struct collSpec collectorData = {
    .numClients = numClients,
    .numServers = numColl,
    .sendRPCData = 1,
  };
  struct clientBufSize bufSizes[numClients_];
  if (tag == STREAMDEFVLIST)
    {
      collectorData.partDesc = NULL;
      collectorData.conversion = NULL;
      bufSizes[0] = computeClientStreamBufSize(serverStreamID, &collectorData);
      collectorData.sendRPCData = 0;
      for (size_t clientIdx = 1; clientIdx < (size_t) numClients_; ++clientIdx)
        bufSizes[clientIdx] = computeClientStreamBufSize(serverStreamID, &collectorData);
    }
  else /* tag == STREAM_DEF_DECOMPOSED_VLIST */
    {
      size_t nVars = (size_t) (vlistNvars(serverVlistID));
      /* unpack client data representation */
      int *conversion = Malloc(nVars * sizeof(*conversion));
      xmpi(MPI_Unpack(buffer, size, pos, conversion, (int) nVars, MPI_INT, pioInterComm));
      collectorData.conversion = rxWin[streamIdx].clientDeco.conversion = conversion;
      /* unpack per-client index lists */
      Xt_idxlist(*partDescPreset)[numClients_] = Malloc(nVars * sizeof(*partDescPreset)),
      *clientPartDesc = Malloc(nVars * sizeof(*clientPartDesc));
      Xt_uid(*partDescUID)[numClients_] = Malloc(nVars * sizeof(*partDescUID));
      rxWin[streamIdx].clientDeco.lists = &partDescPreset[0][0];
      rxWin[streamIdx].clientDeco.uids = &partDescUID[0][0];
      int remainingSize = size - *pos;
      for (int clientIdx = 0; clientIdx < numClients_; ++clientIdx)
        {
          unsigned char *clientBuf = (unsigned char *) buffer + *pos;
          int clientPos = 0;
          for (size_t varIdx = 0; varIdx < nVars; ++varIdx)
            {
              xmpi(MPI_Unpack(clientBuf, remainingSize, &clientPos, partDescUID[varIdx] + clientIdx, 1, YAXT_UID_DT, pioInterComm));
              Xt_uid uid = partDescUID[varIdx][clientIdx];
              for (size_t prevVarIdx = 0; prevVarIdx < varIdx; ++prevVarIdx)
                if (partDescUID[prevVarIdx][clientIdx] == uid)
                  {
                    clientPartDesc[varIdx] = clientPartDesc[prevVarIdx];
                    goto got_idxlist;
                  }
              clientPartDesc[varIdx] = xt_idxlist_unpack(clientBuf, remainingSize, &clientPos, pioInterComm);
            got_idxlist:
              partDescPreset[varIdx][clientIdx] = clientPartDesc[varIdx];
            }
          remainingSize -= clientPos;
          *pos += clientPos;
          collectorData.partDesc = clientPartDesc;
          bufSizes[clientIdx] = computeClientStreamBufSize(serverStreamID, &collectorData);
          collectorData.sendRPCData = 0;
        }
      Free(clientPartDesc);
    }
  createClientStreamBuf(streamIdx, bufSizes, conf);

  size_t maxNumStreamWrites = (size_t) getMaxNumStreamWrites(streamptr),
         currentDstIdxlistCacheSize = cdiPioIdxlistCacheGetSize(DstIdxlistCache);
  if ((neededDstIdxlistCacheSize += maxNumStreamWrites) > currentDstIdxlistCacheSize)
    DstIdxlistCache = cdiPioIdxlistCacheResize(DstIdxlistCache, neededDstIdxlistCacheSize);
  if (conf->cacheXmaps)
    {
      size_t currentXmapCacheSize = cdiPioXmapCacheGetSize(XmapCache);
      if ((neededXmapCacheSize += maxNumStreamWrites) > currentXmapCacheSize)
        XmapCache = cdiPioXmapCacheResize(XmapCache, neededXmapCacheSize);
      if (tag == STREAM_DEF_DECOMPOSED_VLIST)
        buildDecoPresetXmaps(serverStreamID, rxWin[streamIdx].clientDeco, commInqCommColl(), conf);
    }
}

/**
 * @brief is encapsulated in CDI library and run on I/O PEs.
 */

void
cdiPioCollectorMessageLoop()
{
  MPI_Status status;

  xdebug("%s", "START");

  MPI_Comm pioInterComm = cdiPioInqInterComm();
  namespaceSwitchSet(NSSWITCH_STREAM_OPEN_BACKEND, NSSW_FUNC(cdiPioServerStreamOpen));
  namespaceSwitchSet(NSSWITCH_STREAM_CLOSE_BACKEND, NSSW_FUNC(cdiPioServerStreamClose));
  const struct cdiPioConf *conf = cdiPioGetConf();
#ifdef HAVE_PARALLEL_NC4
  cdiPioEnableNetCDFParAccess();
  numPioPrimes = PPM_prime_factorization_32((uint32_t) commInqSizeColl(), &pioPrimes);
#ifdef HAVE_H5GET_LIBVERSION
  {
    unsigned majnum, minnum, relnum;
    extern int H5get_libversion(unsigned *, unsigned *, unsigned *);
    H5get_libversion(&majnum, &minnum, &relnum);
    parH5ZeroCountProblem = (majnum == 1 && ((minnum == 8 && relnum <= 8) || (minnum < 8)));
  }
#endif
#endif
#ifdef HAVE_LIBNETCDF
  cdiSerialOpenFileCount = Calloc(sizeof(cdiSerialOpenFileCount[0]), (size_t) commInqSizeColl());
  namespaceSwitchSet(NSSWITCH_CDF_DEF_TIMESTEP, NSSW_FUNC(cdiPioCdfDefTimestep));
  namespaceSwitchSet(NSSWITCH_CDF_STREAM_SETUP, NSSW_FUNC(cdiPioServerCdfDefVars));
#endif

  int *streamActivity = NULL;
  setupClientRanks();
  if (conf->cacheXmaps) XmapCache = cdiPioXmapCacheNew(16, (size_t) numClients_ + 1);
  for (;;)
    {
      xmpi(MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, pioInterComm, &status));

      int source = status.MPI_SOURCE;
      int tag = status.MPI_TAG;

      switch (tag)
        {
        case FINALIZE:
          xmpi(MPI_Recv(NULL, 0, MPI_INT, source, tag, pioInterComm, &status));
          for (size_t streamIdx = 0; streamIdx < openStreams.size; ++streamIdx)
            {
              int streamID = openStreams.entries[streamIdx];
              if (streamID != CDI_UNDEFID) streamClose(streamID);
            }
          if (rxWin)
            {
              Free(rxWin[0].clientBuf);
              Free(rxWin);
            }
          idSetDestroy(&openStreams);
          Free(aggBuf.mem);
          Free(streamActivity);
          Free(clientRanks_);
#ifdef HAVE_PARALLEL_NC4
          Free(pioPrimes);
#elif defined(HAVE_LIBNETCDF)
          Free(cdiSerialOpenFileCount);
#endif
          if (conf->cacheXmaps) cdiPioXmapCacheDelete(XmapCache);
          cdiPioIdxlistCacheDelete(DstIdxlistCache);
          xdebug("%s", "RETURN");
          return;
        case RESOURCES:
        case STREAMOPEN:
        case STREAMFLUSHCLOSE:
        case STREAMCLOSE:
        case STREAMDEFVLIST:
        case STREAM_DEF_DECOMPOSED_VLIST:
          {
            int size;
            xmpi(MPI_Get_count(&status, MPI_PACKED, &size));
            char *buffer = Malloc((size_t) size);
            xmpi(MPI_Recv(buffer, size, MPI_PACKED, source, tag, pioInterComm, &status));
            int pos = reshUnpackResources(buffer, size, &pioInterComm, (cdiPostResUpdateHook) 0);
            switch (tag)
              {
              case STREAMOPEN: cdiPioRecvStreamOpen(buffer, size, &pos, pioInterComm); break;
              case STREAMFLUSHCLOSE:
              case STREAMCLOSE: cdiPioRecvStreamClose(buffer, size, &pos, pioInterComm, tag == STREAMFLUSHCLOSE, conf); break;
              case STREAMDEFVLIST:
              case STREAM_DEF_DECOMPOSED_VLIST: cdiPioRecvStreamDefVlist(buffer, size, &pos, pioInterComm, tag, conf); break;
              }
            Free(buffer);
            streamActivity = Realloc(streamActivity, openStreams.size * sizeof(streamActivity[0]));
          }
          break;
        case WRITETS:
          {
            xmpi(MPI_Recv(streamActivity, (int) openStreams.size, MPI_INT, source, tag, pioInterComm, &status));
            xdebug("RECEIVED MESSAGE WITH TAG \"WRITETS\": source=%d", source);
            getTimeStepData(streamActivity, conf);
            conf->callbacks[CDIPIO_CALLBACK_POSTWRITEBATCH]();
          }
          break;

        default: xabort("TAG NOT DEFINED!");
        }
    }
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
