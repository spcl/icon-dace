/*
   todo
   README: specialRank Pe closes down, when all output files are closed
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <unistd.h>

#ifdef _POSIX_ASYNCHRONOUS_IO
#include <aio.h>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <mpi.h>

#include "pio.h"
#include "pio_comm.h"
#include "pio_dbuffer.h"
#include "pio_impl.h"
#include "pio_util.h"
#include "dmemory.h"

#define BUF_UNUSED ((unsigned) -1)
struct fileFunnelAIO
{
  char *name;
  struct dBuffer fb;
  struct aiocb *ctrlBlks;
  struct aiocb **cBP;
  /* if bufAssign[i] == BUF_UNUSED, the buffer comprising the bytes in
   * fb.buffer[i*writeAggBufLim..(i+1)*writeAggBufLim) are
   * available, values 0..aioQueueDepth-1 denote active aio_write
   * buffers and higher values in range
   * aioQueueDepth..(aioQueueDepth+numCollectors-1) denote active
   * receive buffers
   */
  unsigned *bufAssign;
  off_t *perCollAmounts, offset, syncOffset;
  unsigned aioQueueDepth;
  unsigned numQueuedWrites;
  int activeCollectors, outstandingSync;
  int handle;
};

/***************************************************************/

static void
initBFiledataPA(struct fileFunnelAIO *bfd, const char *filename, const struct cdiPioConf *conf, size_t nc)
{
  size_t bufSize = conf->writeAggBufLim;
  size_t aioQueueDepth = conf->aioQueueDepth;
  xdebug("filename=%s, buffersize=%zu, ncollectors=%zu, AIO queue depth=%zu", filename, bufSize, nc, aioQueueDepth);

  {
    size_t nameSize = strlen(filename) + 1;
    off_t *restrict perCollAmounts = bfd->perCollAmounts = Malloc(nc * sizeof(*perCollAmounts) * 2 + nameSize);
    char *name = bfd->name = (char *) (bfd->perCollAmounts + nc * 2);
    memcpy(name, filename, nameSize);
    for (size_t i = 0; i < nc; ++i)
      {
        perCollAmounts[i] = 0;
        perCollAmounts[nc + i] = -1;
      }
  }

  if ((bfd->handle = open(bfd->name, O_CREAT | O_WRONLY, 0666)) == -1) xabort("Failed to open %s", bfd->name);

  size_t numBuf = nc + aioQueueDepth;
  cdiPioDbufferInit(&bfd->fb, numBuf * bufSize);
  unsigned *bufAssign = bfd->bufAssign = Malloc(numBuf * sizeof(*bufAssign));
  for (size_t i = 0; i < numBuf; ++i) bufAssign[i] = BUF_UNUSED;
  struct aiocb *ctrlBlks = bfd->ctrlBlks = Calloc(aioQueueDepth, sizeof(*ctrlBlks)),
               **cBP = bfd->cBP = Malloc(aioQueueDepth * sizeof(*cBP));

  for (size_t i = 0; i < aioQueueDepth; i++)
    {
      ctrlBlks[i].aio_fildes = bfd->handle;
      ctrlBlks[i].aio_reqprio = 0;
      ctrlBlks[i].aio_sigevent.sigev_notify = SIGEV_NONE;
      cBP[i] = ctrlBlks + i;
    }

  bfd->aioQueueDepth = (unsigned) aioQueueDepth;
  bfd->numQueuedWrites = 0;
  bfd->offset = 0;
  bfd->syncOffset = 0;
  bfd->activeCollectors = (int) nc;
  bfd->outstandingSync = (int) nc;

  xdebug("filename=%s, opened file, return", bfd->name);
}

/***************************************************************/

static int
destroyBFiledataPA(struct fileFunnelAIO *bfd)
{
  int iret = 0;
  size_t aioQueueDepth = bfd->aioQueueDepth;
  size_t numQueuedWrites = bfd->numQueuedWrites;

  xdebug("filename=%s, cleanup and close file", bfd->name);

  /* close file */
  struct aiocb *ctrlBlks = bfd->ctrlBlks, **cBP = bfd->cBP;
  while (numQueuedWrites > 0)
    {
      size_t numCBP = 0;
      for (size_t i = 0; i < aioQueueDepth; ++i)
        if (ctrlBlks[i].aio_buf != NULL) cBP[numCBP++] = ctrlBlks + i;
      assert(numCBP == numQueuedWrites);
      xdebug("file: %s, numQueuedWrites=%zu", bfd->name, numQueuedWrites);
      do
        {
          iret = aio_suspend((const struct aiocb **) cBP, (int) numQueuedWrites, NULL);
          if (iret < 0 && errno != EINTR) xabort("aio_suspend(2) failed");
        }
      while (iret != 0);
      for (size_t i = 0; i < numQueuedWrites; ++i)
        if ((iret = aio_error(cBP[i])) == 0)
          {
            if (aio_return(cBP[i]) == -1)
              {
                errno = iret;
                SysError("aio_write(2) failed");
              }
            --numQueuedWrites;
            cBP[i]->aio_buf = NULL;
          }
        else if (iret > 0 && iret != EINPROGRESS)
          SysError("aio_write(2) failed");
    }
  bfd->numQueuedWrites = (unsigned) numQueuedWrites;

  if ((iret = ftruncate(bfd->handle, bfd->offset)) == -1) xabort("failed to truncate file %s: %s", bfd->name, strerror(errno));
  if ((iret = close(bfd->handle)) == -1)
    {
      iret = errno;
      Warning("failed to close %s", bfd->name);
    }

  /* file closed, cleanup */

  cdiPioDbufferDestroy(&bfd->fb);
  Free(bfd->perCollAmounts);
  bfd->name = NULL;
  Free(bfd->ctrlBlks);
  Free(bfd->cBP);
  Free(bfd->bufAssign);

  xdebug("%s", "closed file and cleaned up, return");

  return iret;
}

struct awComplete
{
  size_t bufIdx, cbIdx;
};

static struct awComplete
completeSomeAIOWrites(struct fileFunnelAIO *bfd, const struct cdiPioConf *conf)
{
  size_t aioQueueDepth = bfd->aioQueueDepth, numQueuedWrites = bfd->numQueuedWrites;
  assert(aioQueueDepth == numQueuedWrites);
  /* complete aio_write of at least one control block */
  int iret;
  do
    {
      iret = aio_suspend((const struct aiocb **) bfd->cBP, (int) aioQueueDepth, NULL);
      if (iret < 0 && errno != EINTR) xabort("aio_suspend(2) failed");
    }
  while (iret != 0);
  struct aiocb *ctrlBlks = bfd->ctrlBlks;
  size_t bufSize = conf->writeAggBufLim;
  /* only initialized to silence the compiler, the semantics of
   * aio_suspend guarantee for the loop to find at least one
   * completed operation at this point. */
  size_t cbBufIdx = SIZE_MAX, cbIdx = SIZE_MAX;
  /* inspect all control blocks for completion */
  unsigned *bufAssign = bfd->bufAssign;
  for (size_t i = 0; i < aioQueueDepth; ++i)
    if ((iret = aio_error(ctrlBlks + i)) == 0)
      {
        if (aio_return(ctrlBlks + i) == -1)
          {
            errno = iret;
            SysError("failed aio_write");
          }
        --numQueuedWrites;
        cbIdx = i;
        cbBufIdx = (size_t) (ptrdiff_t) ((unsigned char *) ctrlBlks[i].aio_buf - (unsigned char *) bfd->fb.buffer) / bufSize;
        bufAssign[cbBufIdx] = BUF_UNUSED;
        ctrlBlks[i].aio_buf = NULL;
      }
    else if (iret > 0 && iret != EINPROGRESS)
      SysError("failed aio_write");
  bfd->numQueuedWrites = (unsigned) numQueuedWrites;
  return (struct awComplete){ .bufIdx = cbBufIdx, .cbIdx = cbIdx };
}

/***************************************************************/

static size_t
writePA(struct fileFunnelAIO *bfd, int source, size_t nProcsColl, size_t amount, const struct cdiPioConf *conf)
{
  size_t aioQueueDepth = bfd->aioQueueDepth, numQueuedWrites = bfd->numQueuedWrites, numBuf = aioQueueDepth + nProcsColl;

  xdebug("file %s, in", bfd->name);

  struct aiocb *ctrlBlks = bfd->ctrlBlks;
  struct awComplete bufReUse;
  unsigned *bufAssign = bfd->bufAssign;
  if (numQueuedWrites >= aioQueueDepth)
    bufReUse = completeSomeAIOWrites(bfd, conf);
  else
    {
      /* silence compiler uninitialized warning */
      bufReUse = (struct awComplete){ SIZE_MAX, SIZE_MAX };
      for (size_t i = 0; i < aioQueueDepth; ++i)
        if (ctrlBlks[i].aio_buf == NULL)
          {
            bufReUse.cbIdx = i;
            break;
          }
      for (size_t i = 0; i < numBuf; ++i)
        if (bufAssign[i] == BUF_UNUSED)
          {
            bufReUse.bufIdx = i;
            break;
          }
    }

  size_t sourceBufIdx = SIZE_MAX;
  for (size_t i = 0; i < numBuf; ++i)
    if (bufAssign[i] == (unsigned) ((unsigned) source + aioQueueDepth))
      {
        sourceBufIdx = i;
        break;
      }
  assert(sourceBufIdx != SIZE_MAX);
  size_t cbIdx = bufReUse.cbIdx, bufSize = conf->writeAggBufLim;
  ctrlBlks[cbIdx].aio_buf = bfd->fb.buffer + bufSize * sourceBufIdx;
  ctrlBlks[cbIdx].aio_nbytes = amount;
  ctrlBlks[cbIdx].aio_offset = bfd->offset;
  bufAssign[sourceBufIdx] = (unsigned) cbIdx;

  xdebug("before aio_write(), file %s, aio_nbytes=%zu, aio_offset=%zu", bfd->name, ctrlBlks[cbIdx].aio_nbytes,
         ctrlBlks[cbIdx].aio_offset);

  ssize_t iret = aio_write(bfd->ctrlBlks + cbIdx);

  xdebug("after aio_write(), file %s, aio_nbytes=%zu, aio_offset=%zu,"
         "iret=aio_write()=%d",
         bfd->name, bfd->ctrlBlks[cbIdx].aio_nbytes, bfd->ctrlBlks[cbIdx].aio_offset, (int) iret);

  if (iret == -1)
    {
      SysError("did not succeed writing buffer");
    }
  else
    xdebug("buffer written to %s", bfd->name);

  bfd->offset += (off_t) amount;
  bfd->perCollAmounts[source] += (off_t) amount;
  ++bfd->numQueuedWrites;

  xdebug("filename=%s, numQueuedWrites=%u, return", bfd->name, bfd->numQueuedWrites);
  return bufReUse.bufIdx;
}

static size_t
getNextFreeBuf(struct fileFunnelAIO *bfd, size_t nProcsColl, const struct cdiPioConf *conf)
{
  size_t aioQueueDepth = bfd->aioQueueDepth, numQueuedWrites = bfd->numQueuedWrites, numBuf = aioQueueDepth + nProcsColl;
  unsigned *bufAssign = bfd->bufAssign;
  if (numQueuedWrites < aioQueueDepth)
    for (size_t i = 0; i < numBuf; ++i)
      if (bufAssign[i] == BUF_UNUSED) return i;
  return completeSomeAIOWrites(bfd, conf).bufIdx;
}

/***************************************************************/

enum
{
  /* meta operation indices */
  MOOpenRx,
  MOSyncRx,
  MOSyncTx,
  numMetaOp = 3
};

static void
sendCollSync(struct fileFunnelAIO *bfd, int fileID, MPI_Comm commPio, int nProcsColl, MPI_Request (*req)[nProcsColl],
             size_t *numOpenSends, MPI_Status *statui)
{
  if (*numOpenSends) /* make sure previous sync sends completed */
    xmpiStats(MPI_Waitall(nProcsColl, req[MOSyncTx], statui), nProcsColl, statui);
  int sync2CollTag = encodeFileOpTag(fileID, IO_Sync_file);
  for (size_t i = 0; i < (size_t) nProcsColl; ++i)
    xmpi(MPI_Isend(NULL, 0, MPI_INT, (int) i, sync2CollTag, commPio, req[MOSyncTx] + i));
  *numOpenSends = (size_t) (bfd->outstandingSync = nProcsColl);
  bfd->syncOffset = 0;
}

static void
funnelFileCleanup(struct fileFunnelAIO *openFileFunnels, int fileID, size_t *openFileFunnelsFill, size_t *maxOpenFileIDp1,
                  size_t *numOpenRequests, MPI_Comm commPio, int nProcsColl)
{
  xdebug("all are finished with file %d, delete file table entry", fileID);
  struct fileFunnelAIO *bfd = openFileFunnels + fileID;
  assert(bfd->syncOffset == bfd->offset);
  int retval = destroyBFiledataPA(bfd);
  --(*openFileFunnelsFill);
  if ((size_t) fileID + 1 == *maxOpenFileIDp1)
    {
      size_t j = (size_t) fileID + 1;
      while (j > 0 && !openFileFunnels[j - 1].name) --j;
      *maxOpenFileIDp1 = j;
      *numOpenRequests = (size_t) nProcsColl * (numMetaOp + *maxOpenFileIDp1);
    }
  xmpi(MPI_Bcast(&retval, 1, MPI_INT, nProcsColl, commPio));
}

static inline void
reinstallListenReq(struct fileFunnelAIO *bfd, int fileID, int source, size_t nProcsColl, size_t bufSize, unsigned aioQueueDepth,
                   MPI_Request (*req)[nProcsColl], MPI_Comm commPio, const struct cdiPioConf *conf)
{
  size_t bufIdx = getNextFreeBuf(bfd, nProcsColl, conf);
  int tag = encodeFileOpTag(fileID, IO_Send_buffer);
  xmpi(MPI_Irecv(bfd->fb.buffer + bufIdx * bufSize, (int) bufSize, MPI_UNSIGNED_CHAR, source, tag, commPio,
                 req[numMetaOp + fileID] + source));
  bfd->bufAssign[bufIdx] = aioQueueDepth + (unsigned) source;
}

void
pioWriterAIO(void)
{
  size_t openFileFunnelsSize = 1, openFileFunnelsFill = 0;
  struct fileFunnelAIO *openFileFunnels = Malloc(sizeof(*openFileFunnels) * openFileFunnelsSize);
  for (size_t i = 0; i < openFileFunnelsSize; ++i) openFileFunnels[i].name = NULL;
  MPI_Comm commPio = commInqCommPio();
  int nProcsColl = commInqSizeColl();
  const struct cdiPioConf *conf = cdiPioGetConf();
  unsigned aioQueueDepth = conf->aioQueueDepth;
  size_t bufSize = conf->writeAggBufLim;

  assert(aioQueueDepth >= 1);
  xdebug("nProcsColl=%d ", nProcsColl);

  int outstandingFinalizations = nProcsColl;

  int intPackSize;
  xmpi(MPI_Pack_size(2, MPI_INT, commPio, &intPackSize));
  size_t maxPathLen = (size_t) conf->maxPathLen, openMsgSize = (size_t) intPackSize + maxPathLen;
  assert(openMsgSize <= INT_MAX);
  unsigned char *restrict openReqBuffer = Malloc((size_t) nProcsColl * openMsgSize + maxPathLen + 1);
  char *filename = (char *) (openReqBuffer + (size_t) nProcsColl * openMsgSize);
  struct syncMsg *syncMsgs = Malloc((size_t) nProcsColl * sizeof(*syncMsgs));
  /* array of requests first one open request per collector, then one
   * sync receive (close, finalize) request per collector and one for
   * sending sync messages, i.e. 3 rows of fixed meaning followed by
   * one request per file and collector, i.e. number of open files
   * rows.  This array essentially describes all the messages
   * concurrently in flight: open/sync/close/shutdown requests and
   * data record aggregates for each file */
  MPI_Request(*req)[nProcsColl] = Malloc(sizeof(*req) * (numMetaOp + openFileFunnelsSize));
  MPI_Status *statui = Malloc(sizeof(*statui) * (size_t) nProcsColl * (numMetaOp + openFileFunnelsSize)),
             *wstats = Malloc(sizeof(*wstats) * (size_t) nProcsColl);
  int numCompleted, *completed = Malloc(sizeof(*completed) * (size_t) nProcsColl * (numMetaOp + openFileFunnelsSize));
  size_t numOpenSends = 0, maxOpenFileIDp1 = 0, numOpenRequests = (size_t) nProcsColl * (numMetaOp + maxOpenFileIDp1);
  for (size_t i = 0; i < (size_t) nProcsColl; ++i)
    {
      xmpi(MPI_Irecv(openReqBuffer + i * openMsgSize, (int) openMsgSize, MPI_PACKED, (int) i, IO_Open_file, commPio,
                     req[MOOpenRx] + i));
      xmpi(MPI_Irecv(syncMsgs + i, 1, cdiPioSyncMsgDt, (int) i, IO_Sync_file, commPio, req[MOSyncRx] + i));
      req[MOSyncTx][i] = MPI_REQUEST_NULL;
    }
  for (size_t j = 0; j < openFileFunnelsSize; ++j)
    for (size_t i = 0; i < (size_t) nProcsColl; ++i) req[numMetaOp + j][i] = MPI_REQUEST_NULL;

  for (;;)
    {
      assert(numOpenRequests <= INT_MAX);
      xmpiStats(MPI_Waitsome((int) numOpenRequests, *req, &numCompleted, completed, statui), numCompleted, statui);
      for (size_t cmpltIdx = 0; cmpltIdx < (size_t) numCompleted; ++cmpltIdx)
        {
          int rcvIdx = completed[cmpltIdx];
          int source = rcvIdx % nProcsColl;
          int command;
          switch (rcvIdx / nProcsColl)
            {
            case MOOpenRx: command = IO_Open_file; break;
            case MOSyncRx: command = syncMsgs[source].command; break;
            case MOSyncTx: --numOpenSends; continue;
            default: command = IO_Send_buffer;
            }

          xdebug("receive message from source=%d, command=%d (%s)", source, command, cdiPioCmdStrTab[command]);

          switch (command)
            {
            case IO_Open_file:
              {
                int openMsgHdr[2];
                int position = 0;
                unsigned char *reqBuffer = openReqBuffer + (size_t) source * openMsgSize;
                xmpi(MPI_Unpack(reqBuffer, (int) openMsgSize, &position, openMsgHdr, 2, MPI_INT, commPio));
                size_t len = (size_t) openMsgHdr[1];
                assert(len <= maxPathLen);
                xmpi(MPI_Unpack(reqBuffer, (int) openMsgSize, &position, filename, openMsgHdr[1], MPI_CHAR, commPio));
                filename[len] = '\0';
                size_t buffersize = conf->writeAggBufLim;
                xdebug("command %s, filename=%s, buffersize=%zu", cdiPioCmdStrTab[command], filename, buffersize);

                int fileID = openMsgHdr[0];
                size_t prevOpenFileFunnelsSize = openFileFunnelsSize;
                if (openFileFunnelsSize <= (unsigned) fileID)
                  {
                    while (openFileFunnelsSize <= (unsigned) fileID) openFileFunnelsSize *= 2;
                    if (openFileFunnelsSize > (unsigned) INT_MAX + 1) openFileFunnelsSize = (unsigned) INT_MAX + 1;
                    openFileFunnels = Realloc(openFileFunnels, sizeof(*openFileFunnels) * openFileFunnelsSize);
                    for (size_t i = prevOpenFileFunnelsSize; i < openFileFunnelsSize; ++i) openFileFunnels[i].name = NULL;
                    req = Realloc(req, sizeof(*req) * (openFileFunnelsSize + numMetaOp));
                    statui = Realloc(statui, sizeof(*statui) * (size_t) nProcsColl * (openFileFunnelsSize + numMetaOp));
                    completed = Realloc(completed, sizeof(*completed) * (size_t) nProcsColl * (openFileFunnelsSize + numMetaOp));
                    for (size_t j = prevOpenFileFunnelsSize; j < openFileFunnelsSize; ++j)
                      for (size_t i = 0; i < (size_t) nProcsColl; ++i) req[numMetaOp + j][i] = MPI_REQUEST_NULL;
                    maxOpenFileIDp1 = (size_t) fileID + 1;
                  }
                else if (maxOpenFileIDp1 < (size_t) fileID + 1)
                  maxOpenFileIDp1 = (size_t) fileID + 1;
                struct fileFunnelAIO *bfd = openFileFunnels + fileID;
                if (!bfd->name)
                  {
                    for (size_t i = 0; i < prevOpenFileFunnelsSize; ++i)
                      if (openFileFunnels[i].name && !strcmp(openFileFunnels[i].name, filename))
                        xabort("filename %s is already open!"
                               " CDI-PIO does not support concurrent access"
                               " through different filehandles.",
                               filename);
                    initBFiledataPA(bfd, filename, conf, (size_t) nProcsColl);
                    ++openFileFunnelsFill;
                  }
                else if (strcmp(filename, bfd->name) != 0)
                  xabort("filename is not consistent, fileID=%d,"
                         " \"%s\" vs. \"%s\"",
                         fileID, filename, bfd->name);
                reinstallListenReq(bfd, fileID, source, (size_t) nProcsColl, bufSize, aioQueueDepth, req, commPio, conf);
                xmpi(MPI_Irecv(reqBuffer, (int) openMsgSize, MPI_PACKED, source, IO_Open_file, commPio, req[MOOpenRx] + source));
                numOpenRequests = (size_t) nProcsColl * (numMetaOp + maxOpenFileIDp1);
              }
              break;

            case IO_Send_buffer:
              {
                int fileID = rcvIdx / nProcsColl - numMetaOp;
                assert(fileID >= 0 && (size_t) fileID < openFileFunnelsSize && openFileFunnels[fileID].name);
                struct fileFunnelAIO *bfd = openFileFunnels + fileID;

                xdebug("command: %s, fileID=%d, name=%s", cdiPioCmdStrTab[command], fileID, bfd->name);

                int messagesize;
                xmpi(MPI_Get_count(statui + cmpltIdx, MPI_UNSIGNED_CHAR, &messagesize));
                size_t amount = (size_t) messagesize;

                writePA(bfd, source, (size_t) nProcsColl, amount, conf);
                if (!bfd->outstandingSync && bfd->syncOffset == bfd->offset)
                  sendCollSync(bfd, fileID, commPio, nProcsColl, req, &numOpenSends, wstats);
                off_t *restrict perCollAmounts = bfd->perCollAmounts;
                /* end of collector stream not yet reached? */
                if (perCollAmounts[source] != perCollAmounts[nProcsColl + source])
                  reinstallListenReq(bfd, fileID, source, (size_t) nProcsColl, bufSize, aioQueueDepth, req, commPio, conf);
                /* end of collector stream reached */
                else if (!--bfd->activeCollectors)
                  funnelFileCleanup(openFileFunnels, fileID, &openFileFunnelsFill, &maxOpenFileIDp1, &numOpenRequests, commPio,
                                    nProcsColl);
              }
              break;

            case IO_Sync_file:
              {
                int fileID = syncMsgs[source].fileID;
                assert(fileID >= 0 && (size_t) fileID < openFileFunnelsSize && openFileFunnels[fileID].name);
                struct fileFunnelAIO *bfd = openFileFunnels + fileID;
                xdebug("COMMAND %s,  FILE%d, SOURCE%d", cdiPioCmdStrTab[command], fileID, source);
                bfd->syncOffset += syncMsgs[source].amount;
                if (!--bfd->outstandingSync && bfd->syncOffset == bfd->offset)
                  sendCollSync(bfd, fileID, commPio, nProcsColl, req, &numOpenSends, wstats);
                xmpi(MPI_Irecv(syncMsgs + source, 1, cdiPioSyncMsgDt, source, IO_Sync_file, commPio, req[MOSyncRx] + source));
              }
              break;

            case IO_Close_file:
              {
                int fileID = syncMsgs[source].fileID;
                assert(fileID >= 0 && (size_t) fileID < openFileFunnelsSize && openFileFunnels[fileID].name);
                struct fileFunnelAIO *bfd = openFileFunnels + fileID;
                xdebug("COMMAND %s,  FILE%d, SOURCE%d", cdiPioCmdStrTab[command], fileID, source);
                bfd->syncOffset += syncMsgs[source].amount;
                off_t *restrict perCollAmounts = bfd->perCollAmounts;
                perCollAmounts[nProcsColl + source] = syncMsgs[source].amount;
                /* did this source collector already send enough data? */
                if (perCollAmounts[source] == perCollAmounts[nProcsColl + source])
                  {
                    size_t fileReqIdx = numMetaOp + (size_t) fileID;
                    assert(req[fileReqIdx][source] != MPI_REQUEST_NULL);
                    xmpi(MPI_Cancel(req[fileReqIdx] + source));
                    xmpi(MPI_Request_free(req[fileReqIdx] + source));
                    if (!--bfd->activeCollectors)
                      funnelFileCleanup(openFileFunnels, fileID, &openFileFunnelsFill, &maxOpenFileIDp1, &numOpenRequests, commPio,
                                        nProcsColl);
                  }
                xmpi(MPI_Irecv(syncMsgs + source, 1, cdiPioSyncMsgDt, source, IO_Sync_file, commPio, req[MOSyncRx] + source));
              }
              break;

            case IO_Finalize:
              if (req[MOSyncTx][source] != MPI_REQUEST_NULL)
                {
                  xmpiStat(MPI_Wait(req[MOSyncTx] + source, wstats), wstats);
                  --numOpenSends;
                }
              if (req[MOOpenRx][source] != MPI_REQUEST_NULL)
                {
                  xmpi(MPI_Cancel(req[MOOpenRx] + source));
                  xmpi(MPI_Request_free(req[MOOpenRx] + source));
                }
              if (!--outstandingFinalizations)
                {
                  if (openFileFunnelsFill)
                    xabort("some files still not closed.");
                  else
                    {
                      xdebug("%s", "all files are finished, destroy file set,"
                                   " return");
                    }
                  Free(completed);
                  Free(wstats);
                  Free(statui);
                  Free(req);
                  Free(syncMsgs);
                  Free(openReqBuffer);
                  Free(openFileFunnels);
                  return;
                }
              break;
            default: xabort("COMMAND NOT IMPLEMENTED");
            }
        }
    }
}

/***************************************************************/

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
