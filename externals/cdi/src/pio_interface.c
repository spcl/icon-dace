#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>
#include <yaxt.h>

#ifdef HAVE_PPM_CORE
#include <ppm/ppm.h>
#endif

#include "cdi.h"
#include "cdipio.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "namespace.h"
#include "pio.h"
#include "pio_client.h"
#include "pio_conf.h"
#include "pio_dist_grid.h"
#include "pio_id_set.h"
#include "pio_impl.h"
#include "pio_serialize.h"
#include "pio_interface.h"
#include "pio_comm.h"
#include "pio_rpc.h"
#include "pio_server.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "vlist.h"
#include "vlist_var.h"

static struct rdmaWin
{
  size_t size;
  unsigned char *buffer, *head;
  MPI_Win win;
  int dictSize, dictDataUsed, dictRPCUsed;
  bool postSet, refuseFuncCall;
  /** data or meta-data updates on this stream are to be communicated */
  bool pendingUpdate;
  struct partDescPreset clientDeco;
} *txWin = NULL;

static struct idList openStreams;

const char *const funcMap[numRPCFuncs] = {
  "streamDefTimestep",
};

float cdiPIOpartInflate_;

static inline void collWait(size_t streamIdx);

static inline void
collWaitAll()
{
  size_t openStreamsSize = openStreams.size;
  for (size_t streamIdx = 0; streamIdx < openStreamsSize; ++streamIdx)
    if (openStreams.entries[streamIdx] != CDI_UNDEFID) collWait(streamIdx);
}

/****************************************************/

void
memcpyPackFunc(void *dataDesc, void *buf, int size, int *pos, void *context)
{
  struct memCpyDataDesc *p = dataDesc;
  (void) context;
  xassert(size >= *pos && (size_t) (size - *pos) >= p->obj_size);
  memcpy((unsigned char *) buf + *pos, p->obj, p->obj_size);
  *pos += (int) p->obj_size;
}

/****************************************************/

static void
modelWinFlushBuffer(size_t streamIdx)
{
  size_t numStreams = openStreams.size;

  xassert(streamIdx < numStreams && openStreams.entries[streamIdx] != CDI_UNDEFID && txWin != NULL
          && txWin[streamIdx].buffer != NULL);
#if 0
  /* for debugging purposes it might be smart to set the whole buffer
   * to 0 */
  size_t clearSize = txWin[streamIdx].size;
  memset(txWin[streamIdx].buffer, 0, clearSize);
#else
  /* but normally, setting the first header record to zero will
   * suffice (see below) */
#endif
  txWin[streamIdx].head = txWin[streamIdx].buffer + (size_t) txWin[streamIdx].dictSize * sizeof(struct winHeaderEntry);
  txWin[streamIdx].refuseFuncCall = 0;
  txWin[streamIdx].dictDataUsed = 1;
  txWin[streamIdx].dictRPCUsed = 0;
  txWin[streamIdx].pendingUpdate = false;
  txWin[streamIdx].postSet = false;
  *(struct winHeaderEntry *) (void *) txWin[streamIdx].buffer = (struct winHeaderEntry){
    .id = HEADERSIZEMARKER,
    .specific.headerSize = { .numDataEntries = 0, .numRPCEntries = 0 },
  };
}

/************************************************************************/

void
cdiPioClientStreamWinInit(int streamID)
{
  size_t streamIdx = insertID(&openStreams, streamID);
  txWin = Realloc(txWin, openStreams.size * sizeof(txWin[0]));
  txWin[streamIdx].win = MPI_WIN_NULL;
  txWin[streamIdx].buffer = NULL;
  txWin[streamIdx].postSet = false;
  txWin[streamIdx].clientDeco.lists = NULL;
  txWin[streamIdx].clientDeco.uids = NULL;
  txWin[streamIdx].clientDeco.conversion = NULL;
}

void
cdiPioClientStreamWinDestroy(int streamID)
{
  /* no operation on any other window must overlap this call for
   * platforms that block the servers in MPI_Win_complete
   */
  collWaitAll();
  size_t streamIdx = indexOfID(&openStreams, streamID);
  if (txWin[streamIdx].clientDeco.lists)
    {
      size_t nVars = (size_t) (vlistNvars(streamInqVlist(streamID)));
      cdiPioDestroyPartDescPreset(1, nVars, &txWin[streamIdx].clientDeco);
    }
  if (txWin[streamIdx].postSet) xmpi(MPI_Win_wait(txWin[streamIdx].win));
  if (txWin[streamIdx].win != MPI_WIN_NULL) xmpi(MPI_Win_free(&txWin[streamIdx].win));
  if (txWin[streamIdx].buffer) xmpi(MPI_Free_mem(txWin[streamIdx].buffer));
  removeID(&openStreams, streamID);
}

void
cdiPioSetStreamPartDescPreset(int streamID, size_t nVars, const Xt_idxlist *restrict partDesc, const int *restrict conversion)
{
  size_t streamIdx = indexOfID(&openStreams, streamID);
  Xt_idxlist *restrict partDescCopy = Malloc(nVars * sizeof(Xt_idxlist));
  Xt_uid *restrict partDescUID = Malloc(nVars * sizeof(Xt_uid));
  for (size_t i = 0; i < nVars; ++i)
    {
      Xt_uid uid = partDescUID[i] = xt_idxlist_get_uid(partDesc[i]);
      for (size_t j = 0; j < i; ++j)
        if (partDescUID[j] != uid)
          ;
        else
          {
            partDescCopy[i] = partDescCopy[j];
            goto next_index_list;
          }
      partDescCopy[i] = xt_idxlist_copy(partDesc[i]);
    next_index_list:;
    }
  txWin[streamIdx].clientDeco.uids = partDescUID;
  txWin[streamIdx].clientDeco.lists = partDescCopy;
  int *restrict conversionCopy = Malloc(nVars * sizeof(*conversion));
  memcpy(conversionCopy, conversion, nVars * sizeof(*conversion));
  txWin[streamIdx].clientDeco.conversion = conversionCopy;
}

struct partDescPreset
cdiPioGetStreamPartDescPreset(int streamID)
{
  size_t streamIdx = indexOfID(&openStreams, streamID);
  return txWin[streamIdx].clientDeco;
}

void
cdiPioClientStreamWinCreate(int streamID, struct collSpec *cspec)
{
  /* no operation on any other window must overlap this call for
   * platforms that block the servers in MPI_Win_complete
   */
  collWaitAll();
  struct clientBufSize bufSize = computeClientStreamBufSize(streamID, cspec);
  MPI_Info no_locks_info;
  xmpi(MPI_Info_create(&no_locks_info));
  xmpi(MPI_Info_set(no_locks_info, "no_locks", "true"));
  size_t streamIdx = indexOfID(&openStreams, streamID);
  int dictSize = bufSize.numDataRecords + bufSize.numRPCRecords;
  size_t streamBufSize = bufSize.bufSize;
  txWin[streamIdx].dictSize = dictSize;
  txWin[streamIdx].size = streamBufSize;
  xmpi(MPI_Alloc_mem((MPI_Aint) streamBufSize, MPI_INFO_NULL, &txWin[streamIdx].buffer));
  txWin[streamIdx].head = txWin[streamIdx].buffer + (size_t) dictSize * sizeof(struct winHeaderEntry);
  MPI_Comm collClientIntraComm = cdiPioInqCollClientIntraComm();
  xmpi(MPI_Win_create(txWin[streamIdx].buffer, (MPI_Aint) streamBufSize, 1, no_locks_info, collClientIntraComm,
                      &txWin[streamIdx].win));
  modelWinFlushBuffer(streamIdx);
  xmpi(MPI_Info_free(&no_locks_info));
}

bool
cdiPioClientStreamNeedsFlush(int streamID)
{
  size_t streamIdx = indexOfID(&openStreams, streamID);
  xassert(streamIdx != SIZE_MAX);
  return txWin[streamIdx].pendingUpdate & !txWin[streamIdx].postSet;
}

/************************************************************************/

static void
modelWinEnqueue(size_t streamIdx, struct winHeaderEntry header, const void *data, valPackFunc packFunc)
{
  struct winHeaderEntry *winDict = (struct winHeaderEntry *) (void *) txWin[streamIdx].buffer;
  int targetEntry;
  if (header.id > 0 || header.id == PARTDESCMARKER)
    targetEntry = (txWin[streamIdx].dictDataUsed)++;
  else
    targetEntry = txWin[streamIdx].dictSize - ++(txWin[streamIdx].dictRPCUsed);
  if (header.id > 0)
    {
      size_t alignTo = header.id == DATA_HEADER_DOUBLE ? sizeof(double) : sizeof(float);
      int offset = header.offset = (int) roundUpToMultiple((size_t) (txWin[streamIdx].head - txWin[streamIdx].buffer), alignTo);
      MPI_Comm comm = cdiPioInqInterComm();
      packFunc((void *) data, txWin[streamIdx].buffer, (int) txWin[streamIdx].size, &offset, &comm);
      txWin[streamIdx].head = txWin[streamIdx].buffer + offset;
    }
  else if (header.id == PARTDESCMARKER)
    {
      Xt_uid uid = unpackXTUID(header.specific.partDesc.packedUID);
      Xt_idxlist(*restrict partDescPreset) = txWin[streamIdx].clientDeco.lists;
      if (!partDescPreset)
        {
          /* search if same uid entry has already been enqueued */
          for (int entry = 2; entry < targetEntry; entry += 2)
            {
              xassert(winDict[entry].id == PARTDESCMARKER);
              if (unpackXTUID(winDict[entry].specific.partDesc.packedUID) == uid)
                {
                  /* duplicate entries are copied only once per timestep */
                  header.offset = winDict[entry].offset;
                  goto partDescHeaderFinished;
                }
            }
          /* not yet used partition descriptor, serialize at
           * current position */
          int position = header.offset = (int) (txWin[streamIdx].head - txWin[streamIdx].buffer);
          MPI_Comm comm = cdiPioInqInterComm();
          packFunc((void *) data, txWin[streamIdx].buffer, (int) txWin[streamIdx].size, &position, &comm);
          txWin[streamIdx].head = txWin[streamIdx].buffer + position;
        }
      else
        {
          int varID = winDict[targetEntry - 1].specific.dataRecord.varID;
          Xt_uid(*restrict uids) = txWin[streamIdx].clientDeco.uids;
          if (uid != uids[varID] && data != partDescPreset[varID])
            xabort("error: incorrect index list passed to streamWriteVarPart"
                   " or streamWriteVarPartF!");
          header.offset = -1;
        }
    partDescHeaderFinished:;
    }
  else
    {
      int position = header.offset = (int) (txWin[streamIdx].head - txWin[streamIdx].buffer);
      MPI_Comm comm = cdiPioInqInterComm();
      packFunc((void *) data, txWin[streamIdx].buffer, (int) txWin[streamIdx].size, &position, &comm);
      txWin[streamIdx].head = txWin[streamIdx].buffer + position;
    }
  txWin[streamIdx].pendingUpdate = true;
  winDict[targetEntry] = header;
}

static void
cdiPio_xt_idxlist_pack_wrap(void *data, void *buf, int size, int *pos, void *context)
{
  MPI_Comm comm = *(MPI_Comm *) context;
  size_t pack_size = xt_idxlist_get_pack_size((Xt_idxlist) data, comm);
  xassert(size >= *pos && pack_size <= (size_t) (size - *pos));
  xt_idxlist_pack((Xt_idxlist) data, (unsigned char *) buf, size, pos, comm);
}

static inline void
collWait(size_t streamIdx)
{
  if (txWin[streamIdx].postSet)
    {
      xmpi(MPI_Win_wait(txWin[streamIdx].win));
      modelWinFlushBuffer(streamIdx);
    }
}

static inline void
collProbe(size_t streamIdx)
{
  if (!txWin[streamIdx].postSet) return;
  int flag;
  xmpi(MPI_Win_test(txWin[streamIdx].win, &flag));
  if (flag) modelWinFlushBuffer(streamIdx);
}

void
cdiPioRDMAProgress(void)
{
  size_t nStreams = openStreams.size;
  for (size_t streamIdx = 0; streamIdx < nStreams; ++streamIdx)
    if (openStreams.entries[streamIdx] != CDI_UNDEFID) collProbe(streamIdx);
}

static void
cdiPioBufferPartData_(int streamID, int varID, int memtype, const void *packData, valPackFunc packDataFunc, size_t nmiss,
                      Xt_idxlist partDesc)
{
  size_t streamIdx = indexOfID(&openStreams, streamID);
  xassert(streamIdx != SIZE_MAX);
  xassert(varID >= 0 && varID < streamInqNvars(streamID));
  collWaitAll();
  int dataHeaderID = memtype == MEMTYPE_DOUBLE ? DATA_HEADER_DOUBLE : DATA_HEADER_FLOAT;
  struct winHeaderEntry dataHeader = { .id = dataHeaderID, .specific.dataRecord = { varID, nmiss }, .offset = -1 };
  modelWinEnqueue(streamIdx, dataHeader, packData, packDataFunc);
  {
    struct winHeaderEntry partHeader = { .id = PARTDESCMARKER, .offset = 0 };
    Xt_uid uid = xt_idxlist_get_uid(partDesc);
    packXTUID(partHeader.specific.partDesc.packedUID, uid);
    modelWinEnqueue(streamIdx, partHeader, partDesc, cdiPio_xt_idxlist_pack_wrap);
  }
  txWin[streamIdx].refuseFuncCall = 1;
}

static inline size_t
memtype2ElemSize(int memtype)
{
  return memtype == MEMTYPE_DOUBLE ? sizeof(double) : sizeof(float);
}

void
cdiPioBufferPartData(int streamID, int varID, int memtype, const void *data, size_t nmiss, Xt_idxlist partDesc)
{
  size_t chunk = (size_t) (xt_idxlist_get_num_indices(partDesc));
  size_t elemSize = memtype2ElemSize(memtype);
  cdiPioBufferPartData_(streamID, varID, memtype, &(struct memCpyDataDesc){ data, chunk * elemSize }, memcpyPackFunc, nmiss,
                        partDesc);
}

struct scatterGatherDesc
{
  void *data;
  const int *blocklengths, *displacements;
  size_t elemSize;
  unsigned numBlocks;
  unsigned numElems;
};

static void
scatterGatherPackFunc(void *dataDesc, void *buf, int size, int *pos, void *context)
{
  (void) context;
  const struct scatterGatherDesc *p = dataDesc;
  unsigned numBlocks = p->numBlocks;
  const int *bls = p->blocklengths, *disps = p->displacements;
  int pos_ = *pos;
  unsigned char *dstBuf = (unsigned char *) buf + pos_, *bufEnd = (unsigned char *) buf + size;
  size_t elemSize = p->elemSize;
  xassert(elemSize <= SSIZE_MAX);
  const unsigned char *data = p->data;
  unsigned copyCount = 0, numElems = p->numElems;
  for (unsigned j = 0; j < numBlocks && copyCount < numElems; ++j)
    {
      int bl = bls[j];
      if (bl > 0)
        {
          if ((unsigned) bl + copyCount > numElems)
            {
              bl = (int) (numElems - copyCount);
              Warning("%s: %s", "streamWriteScatteredVarPart", "blocks longer than number of elements in index list!");
            }
          size_t bsize = (size_t) bl * elemSize;
          xassert(dstBuf + bsize <= bufEnd);
          memcpy(dstBuf, data + (ssize_t) elemSize * (ssize_t) disps[j], bsize);
          dstBuf += bsize;
        }
    }
  *pos = (int) (dstBuf - (unsigned char *) buf);
}

void
cdiPioBufferPartDataGather(int streamID, int varID, int memtype, const void *data, int numBlocks, const int blocklengths[],
                           const int displacements[], size_t nmiss, Xt_idxlist partDesc)
{
  xassert(numBlocks >= 0);
  cdiPioBufferPartData_(streamID, varID, memtype,
                        &(struct scatterGatherDesc){ .data = (void *) data,
                                                     .blocklengths = blocklengths,
                                                     .displacements = displacements,
                                                     .elemSize = memtype2ElemSize(memtype),
                                                     .numBlocks = (unsigned) numBlocks,
                                                     .numElems = (unsigned) xt_idxlist_get_num_indices(partDesc) },
                        scatterGatherPackFunc, nmiss, partDesc);
}

/************************************************************************/

void
pioBufferFuncCall(int streamID, struct winHeaderEntry header, const void *data, valPackFunc dataPackFunc)
{
  int funcID = header.id;
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);

  xassert(funcID >= MINFUNCID && funcID <= MAXFUNCID);
  xdebug("%s, func: %s", "START", funcMap[(-1 - funcID)]);

  /* FIXME: move this check to cdiPioClientSetup */
  if (clientRank != cdiPioClientRangeStart(collRank, numClients, numColl)) return;

  xassert(txWin != NULL);

  size_t streamIdx = indexOfID(&openStreams, streamID);
  if (streamIdx == SIZE_MAX) streamIdx = insertID(&openStreams, streamID);
  collWaitAll();
  xassert(txWin[streamIdx].dictRPCUsed + txWin[streamIdx].dictDataUsed < txWin[streamIdx].dictSize
          && txWin[streamIdx].refuseFuncCall == 0);
  modelWinEnqueue(streamIdx, header, data, dataPackFunc);
  xdebug("%s", "RETURN");
}

void
cdiPioNoPostCommSetup(void)
{
}

/*****************************************************************************/

static int xtInitByCDI = 0;
#ifdef HAVE_PPM_CORE
static int ppmInitByCDI = 0;
#endif

/* pioInit definition must currently compile even in non-MPI configurations */
/**
   @brief initializes the MPI_Communicators needed for the
  communication between the calculator PEs and the I/O PEs and within the
  group of I/O PEs.

  commGlob: all PEs

  commPIO: I/O PEs, PEs with highest ranks in commGlob

  commModel: calculating PEs, no I/O PEs

  Collective call

  @param commGlob MPI_Communicator of all calling PEs
  @param nProcsIO number of I/O PEs
  @param partInflate allow for array partitions on comute
  PE that are at most sized \f$ partInflate * \lceil arraySize /
  numComputePEs \rceil \f$
  @param postCommSetupActions function which is called by all I/O servers
  after communicator split
  @return int indicating wether the calling PE is a calcutator (1) or not (0)
*/
MPI_Comm
pioInit(MPI_Comm commGlob, int nProcsIO, int IOMode, int *pioNamespace, float partInflate, void (*postCommSetupActions)(void))
{
  int confH = cdiPioConfCreate();
  cdiPioConfSetIOMode(confH, IOMode);
  cdiPioConfSetPartInflate(confH, partInflate);
  cdiPioConfSetCSRole(confH, cdiPioCSRLastN(commGlob, IOMode, nProcsIO));
  cdiPioConfSetPostCommSetupActions(confH, postCommSetupActions);
  MPI_Comm retcomm = cdiPioInit(commGlob, confH, pioNamespace);
  cdiPioConfDestroy(confH);
  return retcomm;
}

static void cdiPioFileWritingInit(void);

MPI_Comm
cdiPioInit(MPI_Comm commGlob, int confResH, int *pioNamespace)
{
  namespaceSwitchSet(NSSWITCH_ABORT, NSSW_FUNC(cdiAbortC_MPI));
  namespaceSwitchSet(NSSWITCH_WARNING, NSSW_FUNC(cdiPioWarning));

  struct cdiPioConf *conf = reshGetVal(confResH, &cdiPioConfOps);

  if (cdiPioExtraNSKeys[cdiPioEKConf] == 0) cdiPioExtraNSKeys[cdiPioEKConf] = cdiNamespaceSwitchNewKey();
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKConf], NSSW_DATA(conf));

  /* check I/O method plausibility */
  int IOMode = conf->IOMode;

  if (IOMode < PIO_MINIOMODE || IOMode > PIO_MAXIOMODE) xabort("IOMODE IS NOT VALID.");

  if (
#ifdef _POSIX_ASYNCHRONOUS_IO
      !sysconf(_SC_ASYNCHRONOUS_IO) &&
#endif
      IOMode == PIO_ASYNCH)
    xabort("CDI-PIO Output method PIO_ASYNCH is unsupported on this system!");

    /* initialize libraries */
#ifdef HAVE_PPM_CORE
  if ((ppmInitByCDI = (!PPM_initialized() || PPM_finalized()))) PPM_initialize(&commGlob, NULL, NULL);
#endif
  if ((xtInitByCDI = (!xt_initialized() || xt_finalized()))) xt_initialize(commGlob);

  int pioNamespace_ = namespaceNew();
  int prevNamespace = namespaceGetActive();
  namespaceSetActive(pioNamespace_);

  /* setup communication */
  cdiPioCommInit(commGlob, IOMode, conf->clientServerRole);
  int sizeGlob = commInqSizeGlob();

  cdiPIOpartInflate_ = conf->partInflate;

  namespaceSetActive(prevNamespace);
  // JUST FOR TEST CASES WITH ONLY ONE MPI TASK
  if (sizeGlob == 1)
    {
      *pioNamespace = pioNamespace_;
      return commInqCommGlob();
    }
  reshRemove(confResH, &cdiPioConfOps);
  namespaceSetActive(pioNamespace_);

#ifdef HAVE_PPM_DIST_ARRAY_H
  reshDistGridUnpack = cdiPioDistGridUnpack;
#endif

  /* check write buffer size consistency */
  {
#ifdef MPI_SIZE_T
    size_t writeAggBufLim = conf->writeAggBufLim;
    xmpi(MPI_Allreduce(MPI_IN_PLACE, &conf->writeAggBufLim, 1, MPI_SIZE_T, MPI_MAX, commGlob));
#else
    unsigned long writeAggBufLim = conf->writeAggBufLim;
    xmpi(MPI_Allreduce(MPI_IN_PLACE, &writeAggBufLim, 1, MPI_UNSIGNED_LONG, MPI_MAX, commGlob));
#endif
    if (conf->writeAggBufLim != writeAggBufLim)
      fprintf(stderr,
              "inconsistent write buffer size detected: %zu vs. "
#ifdef MPI_SIZE_T
              "%lu\n",
#else
              "%zu\n",
#endif
              conf->writeAggBufLim, writeAggBufLim);
#ifndef MPI_SIZE_T
    conf->writeAggBufLim = writeAggBufLim;
#endif
  }
  MPI_Comm modelComm;
  if (commInqIsProcIO())
    {
      int serverNamespace = pioNamespace_;
      namespaceSetActive(serverNamespace);
      namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKConf], NSSW_DATA(conf));
      cdiPioSerializeSetMPI();
      conf->callbacks[CDIPIO_CALLBACK_POSTCOMMSETUP]();
      cdiPioFileWritingInit();
      if (commInqRankColl() >= 0)
        {
          cdiPioCollectorMessageLoop();
          void (*cdiPioFileWritingFinalize)(void) = namespaceSwitchGet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize]).func;
          cdiPioFileWritingFinalize();
        }
      reshListDestruct(serverNamespace);
#ifdef HAVE_PPM_DIST_ARRAY_H
      if (commInqRankColl() >= 0) cdiPioDistGridFinalizeOnce(serverNamespace);
#endif
      cdiPioCommFinalize();
      if (xtInitByCDI) xt_finalize();
#ifdef HAVE_PPM_CORE
      if (ppmInitByCDI) PPM_finalize();
#endif
      modelComm = MPI_COMM_NULL;
      namespaceSetActive(prevNamespace);
      namespaceDelete(serverNamespace);
    }
  else
    {
      *pioNamespace = pioNamespace_;
      cdiPioClientSetup(pioNamespace_, conf);
      modelComm = commInqCommModel();
    }
  namespaceSetActive(prevNamespace);
  reshReplace(confResH, conf, &cdiPioConfOps);
  return modelComm;
}

static void
cdiPioFileWritingFinalizeDefault(void)
{
  xabort("error: failed to setup file writing finalization function!");
}

static void
cdiPioFileWritingInit(void)
{
  if (cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize] == 0)
    cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize] = cdiNamespaceSwitchNewKey();

  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(cdiPioFileWritingFinalizeDefault));
  int IOMode = commInqIOMode();

  xassert(IOMode != PIO_NONE || commInqSizeGlob() == 1);

  switch (IOMode)
    {
    case PIO_NONE: break;
    case PIO_MPI: initMPINONB(); break;
    case PIO_ASYNCH:
    case PIO_WRITER: pioSendInitialize(); break;
    case PIO_FPGUARD: initPOSIXFPGUARDSENDRECV(); break;
    case PIO_MPI_FW_ORDERED: cdiPioFileWriteOrderedInit(); break;
    case PIO_MPI_FW_AT_ALL: cdiPioFileWriteAtAllInit(); break;
    case PIO_MPI_FW_AT_REBLOCK: cdiPioFileWriteAtReblockInit(); break;
    }
}

/*****************************************************************************/

void
pioEndDef(void)
{
  char *buffer;
  int bufferSize;
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);

  if (clientRank == cdiPioClientRangeStart(collRank, numClients, numColl))
    {
      MPI_Comm comm = cdiPioInqInterComm();
      reshPackBufferCreate(&buffer, &bufferSize, &comm);
      xmpi(MPI_Send(buffer, bufferSize, MPI_PACKED, collRank, RESOURCES, comm));
      reshPackBufferDestroy(&buffer);
    }
#if HAVE_PPM_DIST_ARRAY_H
  else
    cdiPioDistGridPackAssist();
#endif
}

/************************************************************************/

void
pioEndTimestepping(void)
{
}

/****************************************************/

/**
 * @brief is invoked by the calculator PEs, to inform
 * the I/O PEs that no more data will be written.
 */
void
pioFinalize(void)
{
  xdebug("%s", "START");

  /* namespace is unchanged on I/O servers */
  int pioNamespace = namespaceGetActive();
  if (pioNamespace == 0) return;
  reshListDestruct(pioNamespace);
#ifdef HAVE_PPM_DIST_ARRAY_H
  cdiPioDistGridFinalizeOnce(pioNamespace);
#endif
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);

  if (clientRank == cdiPioClientRangeStart(collRank, numClients, numColl))
    xmpi(MPI_Send(NULL, 0, MPI_INT, collRank, FINALIZE, cdiPioInqInterComm()));
  xdebug("%s", "SENT MESSAGE WITH TAG \"FINALIZE\"");
  Free(txWin);
  cdiPioCommFinalize();
  idSetDestroy(&openStreams);
  if (xtInitByCDI) xt_finalize();
#ifdef HAVE_PPM_CORE
  if (ppmInitByCDI) PPM_finalize();
#endif
  namespaceDelete(pioNamespace);
  xdebug("%s", "RETURN");
}

/************************************************************************/

static void
cdiPioClientStreamIdxWinPost(size_t streamIdx, MPI_Group remoteGroup)
{
  struct winHeaderEntry header = { .id = HEADERSIZEMARKER,
                                   .specific.headerSize = { .numDataEntries = txWin[streamIdx].dictDataUsed,
                                                            .numRPCEntries = txWin[streamIdx].dictRPCUsed } };
  struct winHeaderEntry *winDict = (struct winHeaderEntry *) (void *) txWin[streamIdx].buffer;
  winDict[0] = header;
  xmpi(MPI_Win_post(remoteGroup, MPI_MODE_NOPUT, txWin[streamIdx].win));
  txWin[streamIdx].postSet = true;
}

void
cdiPioClientStreamWinPost(int streamID)
{
  MPI_Group remoteGroup = cdiPioInqRemoteGroup();
  size_t streamIdx = indexOfID(&openStreams, streamID);
  cdiPioClientStreamIdxWinPost(streamIdx, remoteGroup);
}

void
pioWriteTimestep(void)
{
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);

  xdebug("%s", "START");

  xassert(txWin != NULL);

  size_t nStreams = openStreams.size;
  int *hasUpdates = Malloc(nStreams * sizeof(hasUpdates[0]));
  for (size_t streamIdx = 0; streamIdx < nStreams; ++streamIdx)
    {
      hasUpdates[streamIdx] = txWin[streamIdx].pendingUpdate;
    }

  if (clientRank == cdiPioClientRangeStart(collRank, numClients, numColl))
    {
      xmpi(MPI_Send(hasUpdates, (int) nStreams, MPI_INT, collRank, WRITETS, cdiPioInqInterComm()));
      xdebug("%s", "SENT MESSAGE WITH TAG \"WRITETS\"");
    }

  MPI_Group remoteGroup = cdiPioInqRemoteGroup();
  for (size_t streamIdx = 0; streamIdx < nStreams; ++streamIdx)
    if (openStreams.entries[streamIdx] != CDI_UNDEFID && hasUpdates[streamIdx])
      cdiPioClientStreamIdxWinPost(streamIdx, remoteGroup);

  Free(hasUpdates);
  xdebug("%s", "RETURN. messages sent, windows posted");
}

void
cdiPioStreamWriteVarPart_(int streamID, int varID, int memtype, const void *data, int nmiss, Xt_idxlist partDesc)
{
  if (CDI_Debug) Message("streamID = %d varID = %d", streamID, varID);

  int chunk = xt_idxlist_get_num_indices(partDesc);
  xassert(chunk == 0 || data);

  void (*myStreamWriteVarPart)(int streamID, int varID, int memtype, const void *data, int nmiss, Xt_idxlist partDesc)
      = (void (*)(int, int, int, const void *, int, Xt_idxlist)) namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_PART_).func;

  if (!myStreamWriteVarPart) xabort("local part writing is unsupported!");

  myStreamWriteVarPart(streamID, varID, memtype, data, nmiss, partDesc);
}

void
streamWriteVarPart(int streamID, int varID, const double *data, int nmiss, Xt_idxlist partDesc)
{
  cdiPioStreamWriteVarPart_(streamID, varID, MEMTYPE_DOUBLE, data, nmiss, partDesc);
}

void
streamWriteVarPartF(int streamID, int varID, const float *data, int nmiss, Xt_idxlist partDesc)
{
  cdiPioStreamWriteVarPart_(streamID, varID, MEMTYPE_FLOAT, data, nmiss, partDesc);
}

static void
streamWriteScatteredVarPart_(int streamID, int varID, int memtype, const void *data, int numBlocks, const int blocklengths[],
                             const int displacements[], int nmiss, Xt_idxlist partDesc)
{
  if (CDI_Debug) Message("streamID = %d varID = %d", streamID, varID);

  int chunk = xt_idxlist_get_num_indices(partDesc);
  xassert(chunk == 0 || data);

  void (*myStreamWriteScatteredVarPart)(int streamID, int varID, int memtype, const void *data, int numBlocks,
                                        const int blocklengths[], const int displacements[], int nmiss, Xt_idxlist partDesc)
      = (void (*)(int, int, int, const void *, int, const int[], const int[], int, Xt_idxlist)) namespaceSwitchGet(
            NSSWITCH_STREAM_WRITE_SCATTERED_VAR_PART_)
            .func;

  if (!myStreamWriteScatteredVarPart) xabort("local part writing is unsupported!");

  myStreamWriteScatteredVarPart(streamID, varID, memtype, data, numBlocks, blocklengths, displacements, nmiss, partDesc);
}

void
streamWriteScatteredVarPart(int streamID, int varID, const double *data, int numBlocks, const int blocklengths[],
                            const int displacements[], int nmiss, Xt_idxlist partDesc)
{
  streamWriteScatteredVarPart_(streamID, varID, MEMTYPE_DOUBLE, data, numBlocks, blocklengths, displacements, nmiss, partDesc);
}

void
streamWriteScatteredVarPartF(int streamID, int varID, const float *data, int numBlocks, const int blocklengths[],
                             const int displacements[], int nmiss, Xt_idxlist partDesc)
{
  streamWriteScatteredVarPart_(streamID, varID, MEMTYPE_FLOAT, data, numBlocks, blocklengths, displacements, nmiss, partDesc);
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
