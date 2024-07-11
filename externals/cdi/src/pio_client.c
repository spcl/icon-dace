#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>

#include <mpi.h>
#include <yaxt.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#ifdef HAVE_LIBGRIB
#include "gribapi.h"
#endif
#include "namespace.h"
#include "taxis.h"
#include "vlist.h"

#include "cdipio.h"
#include "pio.h"
#include "pio_client.h"
#include "pio_dist_grid.h"
#include "pio_comm.h"
#include "pio_id_set.h"
#include "pio_idxlist_cache.h"
#include "pio_interface.h"
#include "pio_rpc.h"
#include "pio_util.h"
#include "pio_serialize.h"

struct ResHListUpdate
{
  unsigned char *msgBuffer;
  MPI_Comm comm;
  int numClients, clientRank, numColl, collRank, groupLeader;
  int msgSize, msgPos;
  bool sendRPCData;
};

static struct ResHListUpdate
cdiPioClientUpdateResHList()
{
  MPI_Comm comm = cdiPioInqInterComm();
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);
  int groupLeader = cdiPioClientRangeStart(collRank, numClients, numColl);
  int sendRPCData = (clientRank == groupLeader);
  char *msgBuffer = NULL;
  int msgSize = 0, msgBufPos = 0;
  if (sendRPCData)
    {
      msgBufPos = reshPackBufferCreate(&msgBuffer, &msgSize, &comm);
    }
#ifdef HAVE_PPM_DIST_ARRAY_H
  else
    cdiPioDistGridPackAssist();
#endif
  return (struct ResHListUpdate){ .msgBuffer = (unsigned char *) msgBuffer,
                                  .comm = comm,
                                  .numClients = numClients,
                                  .clientRank = clientRank,
                                  .groupLeader = groupLeader,
                                  .numColl = numColl,
                                  .collRank = collRank,
                                  .msgSize = msgSize,
                                  .msgPos = msgBufPos,
                                  .sendRPCData = sendRPCData };
}

static int
cdiPioClientStreamOpen(const char *filename, char filemode, int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
{
  (void) streamptr;
  (void) recordBufIsToBeCreated;
  int fileID = 0;
  if (filemode == 'w')
    {
      int streamID = streamptr->self;
      reshSetStatus(streamID, &streamOps, reshGetStatus(streamID, &streamOps) & ~RESH_SYNC_BIT);
      struct ResHListUpdate rup = cdiPioClientUpdateResHList();
      MPI_Comm comm = rup.comm;
      int collRank = rup.collRank;
      if (rup.sendRPCData)
        {
          unsigned char *msgBuffer = rup.msgBuffer;
          int msgSize = rup.msgSize;
          int size;
          size_t filename_len = strlen(filename);
          xassert(filename_len < (size_t) (INT_MAX - rup.msgPos));
          int soHdr[3] = { streamptr->self, filetype, (int) filename_len };
          xmpi(MPI_Pack_size(3, MPI_INT, comm, &size));
          msgSize += size;
          xmpi(MPI_Pack_size(1, MPI_CHAR, comm, &size));
          msgSize += size;
          xmpi(MPI_Pack_size((int) filename_len, MPI_CHAR, comm, &size));
          msgSize += size;
          /* optimize to pos + size */
          msgBuffer = Realloc(msgBuffer, (size_t) msgSize);
          xmpi(MPI_Pack(soHdr, 3, MPI_INT, msgBuffer, msgSize, &rup.msgPos, comm));
          xmpi(MPI_Pack(&filemode, 1, MPI_CHAR, msgBuffer, msgSize, &rup.msgPos, comm));
          xmpi(MPI_Pack((void *) filename, (int) filename_len, MPI_CHAR, msgBuffer, msgSize, &rup.msgPos, comm));
          xmpi(MPI_Sendrecv(msgBuffer, msgSize, MPI_PACKED, collRank, STREAMOPEN, &fileID, 1, MPI_INT, collRank, STREAMOPEN, comm,
                            MPI_STATUS_IGNORE));
          Free(msgBuffer);
        }
      else
        {
          xmpi(MPI_Recv(&fileID, 1, MPI_INT, collRank, STREAMOPEN, comm, MPI_STATUS_IGNORE));
        }
      if (fileID >= 0)
        {
          streamptr->filetype = filetype;
          cdiPioClientStreamWinInit(streamID);
        }
    }
  else
    xabort("cdiPIO read support not implemented");
  return fileID;
}

static void
checkVlistForPIO(int vlistID)
{
  int nVars = vlistNvars(vlistID);
  int varShape[3];
  for (int varID = 0; varID < nVars; ++varID)
    {
      cdiPioQueryVarDims(varShape, vlistID, varID);
      intmax_t size = 1;
      for (size_t i = 0; i < 3; ++i) size *= varShape[i];
      if (size > XT_INT_MAX)
        {
          errno = EOVERFLOW;
          SysError("CDI Variable too large for YAXT,"
                   " use configure with larger index type\n"
                   "vlistID=%d, varID=%d, variable size=%zu",
                   vlistID, varID, (size_t) size);
        }
    }
}

static void
cdiPioClientStreamDefVlistCommon(int streamID, int vlistID)
{
  checkVlistForPIO(vlistID);
  cdiStreamDefVlist_(streamID, vlistID);
  reshSetStatus(streamID, &streamOps, reshGetStatus(streamID, &streamOps) & ~RESH_SYNC_BIT);
}

static void
cdiPioClientStreamDefVlist_(int streamID, int vlistID)
{
  cdiPioClientStreamDefVlistCommon(streamID, vlistID);
  struct ResHListUpdate rup = cdiPioClientUpdateResHList();
  if (rup.sendRPCData)
    {
      MPI_Comm comm = rup.comm;
      unsigned char *msgBuffer = rup.msgBuffer;
      /* optimize: pos + size */
      int msgSize = rup.msgSize;
      {
        int size;
        xmpi(MPI_Pack_size(defVlistNInts, MPI_INT, comm, &size));
        msgSize += size;
      }
      msgBuffer = Realloc(msgBuffer, (size_t) msgSize);
      int msgData[defVlistNInts] = { streamID, streamInqVlist(streamID) };
      xmpi(MPI_Pack(&msgData, defVlistNInts, MPI_INT, msgBuffer, msgSize, &rup.msgPos, comm));
      xmpi(MPI_Send(msgBuffer, rup.msgPos, MPI_PACKED, rup.collRank, STREAMDEFVLIST, comm));
      Free(msgBuffer);
    }
  struct collSpec cspec = {
    .numClients = rup.numClients, .numServers = rup.numColl, .sendRPCData = rup.sendRPCData, .partDesc = NULL, .conversion = NULL
  };
  cdiPioClientStreamWinCreate(streamID, &cspec);
}

static int
sizePartDescPresetPack(size_t nVars, const struct partDescPreset clientDeco, bool *isDuplicate, MPI_Comm comm)
{
  int packSize = 0;
  const Xt_uid *uids = clientDeco.uids;
  Xt_idxlist *partDesc = clientDeco.lists;
  {
    int sizeXtUID;
    xmpi(MPI_Pack_size(1, YAXT_UID_DT, comm, &sizeXtUID));
    packSize += (int) nVars * sizeXtUID;
  }
  for (size_t varIdx = 0; varIdx < nVars; ++varIdx)
    {
      Xt_uid uid = uids[varIdx];
      for (size_t j = 0; j < varIdx; ++j)
        if (uids[j] != uid)
          ;
        else
          {
            isDuplicate[varIdx] = true;
            goto next_index_list;
          }
      isDuplicate[varIdx] = false;
      packSize += (int) xt_idxlist_get_pack_size(partDesc[varIdx], comm);
    next_index_list:;
    }
  return packSize;
}

static int
packPartDescPreset(void *packBuf, int bufSize, size_t nVars, struct partDescPreset clientDeco, const bool *isDuplicate,
                   MPI_Comm comm)
{
  const Xt_uid *uids = clientDeco.uids;
  Xt_idxlist *partDesc = clientDeco.lists;
  int position = 0;
  for (size_t varIdx = 0; varIdx < nVars; ++varIdx)
    {
      xmpi(MPI_Pack((void *) (uids + varIdx), 1, YAXT_UID_DT, packBuf, bufSize, &position, comm));
      if (!isDuplicate[varIdx]) xt_idxlist_pack(partDesc[varIdx], packBuf, bufSize, &position, comm);
    }
  return position;
}

void
cdiPioStreamDefDecomposedVlist(int streamID, int vlistID, const Xt_idxlist partDesc[], const int conversion[])
{
  cdiPioClientStreamDefVlistCommon(streamID, vlistID);
  struct ResHListUpdate rup = cdiPioClientUpdateResHList();
  MPI_Comm comm = rup.comm;
  int collRank = rup.collRank;
  size_t nVars = (size_t) (vlistNvars(vlistID));
  MPI_Comm aggComm = cdiPioInqCollClientIntraComm();
  cdiPioAssertConsistentIntVec(nVars, conversion, aggComm, rup.clientRank, 0, rup.numClients);
  cdiPioSetStreamPartDescPreset(streamID, nVars, partDesc, conversion);
  int lastClientInGroup = cdiPioClientRangeStart(collRank + 1, rup.numClients, rup.numColl) - 1, clientRank = rup.clientRank,
      numClientsInGroup = lastClientInGroup - clientRank + 1;
  struct partDescPreset clientDeco = cdiPioGetStreamPartDescPreset(streamID);
  bool *partDescDuplicates = Malloc(nVars * sizeof(*partDescDuplicates));
  int listPackSizeTotal = sizePartDescPresetPack(nVars, clientDeco, partDescDuplicates, comm);
  if (rup.sendRPCData)
    {
      unsigned char *msgBuffer = rup.msgBuffer;
      /* optimize: pos + size */
      int msgSize = rup.msgSize;
      /* TODO: find size of idxlists for packing */
      int(*listPackSizes) = Malloc((size_t) numClientsInGroup * sizeof(*listPackSizes));
      listPackSizes[0] = listPackSizeTotal;
      MPI_Request *req = (numClientsInGroup - 1) ? Malloc((size_t) (numClientsInGroup - 1) * sizeof(*req)) : NULL;
      for (int rank = rup.clientRank + 1; rank <= lastClientInGroup; ++rank)
        {
          size_t rankOfs = (size_t) (rank - clientRank);
          xmpi(MPI_Irecv(listPackSizes + rankOfs, 1, MPI_INT, rank, IDXLIST_PRESET_SIZE, aggComm, req + rankOfs - 1));
        }
      xmpi(MPI_Waitall(numClientsInGroup - 1, req, MPI_STATUSES_IGNORE));
      for (int rank = rup.clientRank; rank <= lastClientInGroup; ++rank)
        {
          size_t rankOfs = (size_t) (rank - clientRank);
          msgSize += listPackSizes[rankOfs];
        }
      {
        int size;
        xmpi(MPI_Pack_size(defVlistNInts + (int) nVars, MPI_INT, comm, &size));
        msgSize += size;
      }
      /* size of extra packaging:  */
      msgBuffer = Realloc(msgBuffer, (size_t) msgSize);
      int msgData[defVlistNInts] = { streamID, streamInqVlist(streamID) };
      xmpi(MPI_Pack(&msgData, defVlistNInts, MPI_INT, msgBuffer, msgSize, &rup.msgPos, comm));
      xmpi(MPI_Pack((int *) conversion, (int) nVars, MPI_INT, msgBuffer, msgSize, &rup.msgPos, comm));
      rup.msgPos += packPartDescPreset(msgBuffer + rup.msgPos, listPackSizeTotal, nVars, clientDeco, partDescDuplicates, comm);
      for (int rank = rup.clientRank + 1; rank <= lastClientInGroup; ++rank)
        {
          size_t rankOfs = (size_t) (rank - clientRank);
          int clientPackSize = listPackSizes[rankOfs];
          xmpi(MPI_Irecv(msgBuffer + rup.msgPos, clientPackSize, MPI_PACKED, rank, IDXLIST_PRESET, aggComm, req + rankOfs - 1));
          rup.msgPos += clientPackSize;
        }
      xmpi(MPI_Waitall(numClientsInGroup - 1, req, MPI_STATUSES_IGNORE));
      xmpi(MPI_Send(msgBuffer, rup.msgPos, MPI_PACKED, collRank, STREAM_DEF_DECOMPOSED_VLIST, comm));
      Free(msgBuffer);
      Free(req);
      Free(listPackSizes);
    }
  else
    {
      int groupLeader = rup.groupLeader;
      MPI_Request req[2];
      unsigned char *ilPackBuf = Malloc((size_t) listPackSizeTotal);
      int position = packPartDescPreset(ilPackBuf, listPackSizeTotal, nVars, clientDeco, partDescDuplicates, comm);
      xmpi(MPI_Isend(&position, 1, MPI_INT, groupLeader, IDXLIST_PRESET_SIZE, aggComm, req));
      xmpi(MPI_Isend(ilPackBuf, position, MPI_PACKED, groupLeader, IDXLIST_PRESET, aggComm, req + 1));
      xmpi(MPI_Waitall(2, req, MPI_STATUSES_IGNORE));
      Free(ilPackBuf);
    }
  Free(partDescDuplicates);
  struct collSpec cspec = { .numClients = rup.numClients,
                            .numServers = rup.numColl,
                            .sendRPCData = rup.sendRPCData,
                            .partDesc = partDesc,
                            .conversion = (int *) conversion };
  cdiPioClientStreamWinCreate(streamID, &cspec);
}

static void
cdiPioClientStreamWriteVar_(int streamID, int varID, int memtype, const void *data, size_t numMissVals)
{
  struct partDescPreset clientDeco = cdiPioGetStreamPartDescPreset(streamID);
  int vlistID = streamInqVlist(streamID);
  size_t nVars = (size_t) (vlistNvars(vlistID));
  if (varID < 0 || (size_t) varID > nVars) xabort("invalid variable ID %d for stream %d requested!", varID, streamID);
  if (!clientDeco.uids)
    xabort("parallel writing requires explicit partition information,"
           " use either streamWriteVarPart or cdiPioStreamDefDecomposedVlist!");
  if (clientDeco.conversion && memtype != clientDeco.conversion[varID])
    xabort("data type does not match pre-declared conversion for stream %d,"
           " variable ID %d!",
           streamID, varID);
  cdiPioStreamWriteVarPart_(streamID, varID, memtype, data, numMissVals, clientDeco.lists[varID]);
}

static struct cdiPioIdxlistCache *clientIdxlistCache;
static size_t neededClientIdxlistCacheSize;
static struct idList seenStreamWriteVarChunk;

static void
cdiPioClientStreamWriteVarChunk_(int streamID, int varID, int memtype, const int rect[][2], const void *data, size_t numMissVals)
{
  if (memtype != MEMTYPE_DOUBLE && memtype != MEMTYPE_FLOAT) Error("Writing of type data %d not implemented!", memtype);
  int vlistID = streamInqVlist(streamID);
  if (indexOfID(&seenStreamWriteVarChunk, vlistID) == SIZE_MAX)
    {
      insertID(&seenStreamWriteVarChunk, vlistID);
      int nvars = vlistNvars(vlistID);
      neededClientIdxlistCacheSize += (size_t) nvars;
      clientIdxlistCache = cdiPioIdxlistCacheResize(clientIdxlistCache, neededClientIdxlistCacheSize);
    }
  int size = vlistInqVarSize(vlistID, varID), varShape[3];
  int ndims = cdiPioQueryVarDims(varShape, vlistID, varID);
  Xt_int varShapeXt[3], origin[3] = { 0, 0, 0 };
  int chunkShape[3] = { 1, 1, 1 };
  if (ndims == 3)
    for (int i = 0; i < 3; ++i) varShapeXt[i] = varShape[i];
  else
    {
      varShapeXt[0] = varShape[0];
      varShapeXt[1] = varShape[2];
    }
  for (int i = 0; i < ndims; ++i) chunkShape[i] = rect[i][1] - rect[i][0] + 1;
  size_t varSize = (size_t) varShape[0] * (size_t) varShape[1] * (size_t) varShape[2];
  xassert(varSize == (size_t) size);
  if (clientIdxlistCache)
    ;
  else
    clientIdxlistCache = cdiPioIdxlistCacheNew(31);
  Xt_idxlist (*cacheSection)(struct cdiPioIdxlistCache * cache, const Xt_int wholeShape[], const Xt_int sliceOrigin[],
                             const int sliceShape[])
      = (ndims == 3) ? cdiPioIdxlistCacheAddSection3D : cdiPioIdxlistCacheAddSection2D;
  Xt_idxlist chunkDesc = cacheSection(clientIdxlistCache, varShapeXt, origin, chunkShape);
  cdiPioBufferPartData(streamID, varID, memtype, data, numMissVals, chunkDesc);
}

#if defined HAVE_LIBNETCDF
static void
cdiPioCdfDefTimestepNOP(stream_t *streamptr, int tsID)
{
  (void) streamptr;
  (void) tsID;
}
#endif

static void
cdiPioClientStreamNOP(stream_t *streamptr)
{
  (void) streamptr;
}

static void
cdiPioClientStreamClose(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  int streamID = streamptr->self;
  int clientRank = commInqRankModel(), numClients = cdiPioCommInqSizeClients(), numColl = commInqSizeColl(),
      collRank = cdiPioCollRank(clientRank, numClients, numColl);
  reshSetStatus(streamID, &streamOps, reshGetStatus(streamID, &streamOps) & ~RESH_SYNC_BIT);
  struct ResHListUpdate rup = cdiPioClientUpdateResHList();
  bool needsFlush = cdiPioClientStreamNeedsFlush(streamID);
  if (rup.sendRPCData)
    {
      MPI_Comm comm = rup.comm;
      unsigned char *msgBuffer = rup.msgBuffer;
      int msgSize = rup.msgSize;
      {
        int size;
        xmpi(MPI_Pack_size(1, MPI_INT, comm, &size));
        msgSize += size;
      }
      /* optimize: pos + size */
      msgBuffer = Realloc(msgBuffer, (size_t) msgSize);
      xmpi(MPI_Pack(&streamptr->self, 1, MPI_INT, msgBuffer, msgSize, &rup.msgPos, comm));
      int tag = needsFlush ? STREAMFLUSHCLOSE : STREAMCLOSE;
      xmpi(MPI_Send(msgBuffer, rup.msgPos, MPI_PACKED, collRank, tag, comm));
      Free(msgBuffer);
    }
  if (needsFlush) cdiPioClientStreamWinPost(streamID);
  cdiPioClientStreamWinDestroy(streamID);
  if (recordBufIsToBeDeleted) switch (streamptr->filetype)
      {
#ifdef HAVE_LIBGRIB
      case CDI_FILETYPE_GRB:
      case CDI_FILETYPE_GRB2: gribContainersDelete(streamptr); break;
#endif
      }
}

static void
cdiPioTaxisPackWrap(void *data, void *buf, int size, int *pos, void *context)
{
  int taxisID = (int) (intptr_t) data;
  reshPackResource(taxisID, &taxisOps, buf, size, pos, context);
}

static int
cdiPioClientStreamDefTimestep_(stream_t *streamptr, int tsID)
{
  int taxisID = vlistInqTaxis(streamptr->vlistID);
  struct winHeaderEntry header = { .id = STREAMDEFTIMESTEP, .specific.funcArgs.streamNewTimestep = { streamptr->self, tsID } };
  xassert(sizeof(void *) >= sizeof(int));
  pioBufferFuncCall(streamptr->self, header, (void *) (intptr_t) taxisID, cdiPioTaxisPackWrap);
  return cdiStreamDefTimestep_(streamptr, tsID);
}

static void
cdiPioClientVlistDestroy_(int vlistID, bool assertInternal)
{
  if (indexOfID(&seenStreamWriteVarChunk, vlistID) != SIZE_MAX)
    {
      int nvars = vlistNvars(vlistID);
      neededClientIdxlistCacheSize -= (size_t) nvars;
      clientIdxlistCache = cdiPioIdxlistCacheResize(clientIdxlistCache, neededClientIdxlistCacheSize);
      removeID(&seenStreamWriteVarChunk, vlistID);
    }
  cdiVlistDestroy_(vlistID, assertInternal);
}

void
cdiPioClientSetup(int pioNamespace, struct cdiPioConf *conf)
{
  int callerCDINamespace = namespaceGetActive();
  namespaceSetActive(pioNamespace);
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKConf], NSSW_DATA(conf));
  cdiPioSerializeSetMPI();
  namespaceSwitchSet(NSSWITCH_ABORT, NSSW_FUNC(cdiAbortC_MPI));
  namespaceSwitchSet(NSSWITCH_WARNING, NSSW_FUNC(cdiPioWarning));
  namespaceSwitchSet(NSSWITCH_STREAM_OPEN_BACKEND, NSSW_FUNC(cdiPioClientStreamOpen));
  namespaceSwitchSet(NSSWITCH_STREAM_DEF_VLIST_, NSSW_FUNC(cdiPioClientStreamDefVlist_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_, NSSW_FUNC(cdiPioClientStreamWriteVar_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_, NSSW_FUNC(cdiPioClientStreamWriteVarChunk_));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_VAR_PART_, NSSW_FUNC(cdiPioBufferPartData));
  namespaceSwitchSet(NSSWITCH_STREAM_WRITE_SCATTERED_VAR_PART_, NSSW_FUNC(cdiPioBufferPartDataGather));
  namespaceSwitchSet(NSSWITCH_STREAM_CLOSE_BACKEND, NSSW_FUNC(cdiPioClientStreamClose));
  namespaceSwitchSet(NSSWITCH_STREAM_DEF_TIMESTEP_, NSSW_FUNC(cdiPioClientStreamDefTimestep_));
  namespaceSwitchSet(NSSWITCH_STREAM_SYNC, NSSW_FUNC(cdiPioClientStreamNOP));
  namespaceSwitchSet(NSSWITCH_VLIST_DESTROY_, NSSW_FUNC(cdiPioClientVlistDestroy_));
#ifdef HAVE_LIBNETCDF
  namespaceSwitchSet(NSSWITCH_CDF_DEF_TIMESTEP, NSSW_FUNC(cdiPioCdfDefTimestepNOP));
  namespaceSwitchSet(NSSWITCH_CDF_STREAM_SETUP, NSSW_FUNC(cdiPioClientStreamNOP));
#endif
  namespaceSetActive(callerCDINamespace);
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
