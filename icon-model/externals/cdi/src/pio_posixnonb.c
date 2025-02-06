#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include "cdi.h"
#include "dmemory.h"

#include "pio.h"
#include "pio_comm.h"
#include "pio_conf.h"
#include "pio_dbuffer.h"
#include "pio_impl.h"
#include "pio_util.h"

struct fileFunnelStdio
{
  char *name;
  struct dBuffer *db;
  FILE *fp;
  off_t *perCollAmounts, offset, syncOffset;
  int activeCollectors, outstandingSync;
};

/***************************************************************/

static void
initBFiledataP(struct fileFunnelStdio *bfp, const char *filename, size_t bs, size_t nc)
{
  xdebug("filename=%s, buffersize=%lu, ncollectors=%zu", filename, bs, nc);

  {
    size_t nameSize = strlen(filename) + 1;
    off_t *restrict perCollAmounts = bfp->perCollAmounts = Malloc(nc * sizeof(*perCollAmounts) * 2 + nameSize);
    char *name = bfp->name = (char *) (bfp->perCollAmounts + nc * 2);
    memcpy(name, filename, nameSize);
    for (size_t i = 0; i < nc; ++i)
      {
        perCollAmounts[i] = 0;
        perCollAmounts[nc + i] = -1;
      }
  }

  if ((bfp->fp = fopen(filename, "w")) == NULL)
    {
      int fopen_errno = errno;
      xabort("Failed to open %s\nerrno=%d: %s", bfp->name, fopen_errno, strerror(fopen_errno));
    }

  struct dBuffer *db = bfp->db = Malloc(sizeof(*db) * nc);
  for (size_t i = 0; i < nc; ++i) cdiPioDbufferInit(db + i, bs);
  bfp->offset = 0;
  bfp->syncOffset = 0;
  bfp->activeCollectors = (int) nc;
  bfp->outstandingSync = (int) nc;
  xdebug("filename=%s, opened file, return", bfp->name);
}

/***************************************************************/

static int
destroyBFiledataP(struct fileFunnelStdio *bfp)
{

  xdebug("filename=%s, cleanup, in", bfp->name);

  /* close file */
  int iret;
  if ((iret = fclose(bfp->fp)) == EOF)
    {
      iret = errno;
      Warning("Failed to close %s", bfp->name);
    }
  /* file closed, cleanup */
  size_t numColl = (size_t) (commInqSizeColl());
  struct dBuffer *db = bfp->db;
  for (size_t i = 0; i < numColl; ++i) cdiPioDbufferDestroy(db + i);
  Free(db);
  Free(bfp->perCollAmounts);
  bfp->name = NULL;

  xdebug("%s", "cleaned up, return");

  return iret;
}

/***************************************************************/

static void
writeP(struct fileFunnelStdio *bfd, int source, size_t amount)
{
  size_t written;
  xdebug("filename=%s, amount=%ld, in", bfd->name, amount);

  bfd->perCollAmounts[source] += (off_t) amount;
  bfd->offset += (off_t) amount;
  if ((written = fwrite(bfd->db[source].buffer, 1, amount, bfd->fp)) != amount)
    xabort("did not succeed writing buffer to %s", bfd->name);

  xdebug("filename=%s, written=%ld, amount=%ld, return", bfd->name, written, amount);
}

enum
{
  /* meta operation indices */
  MOOpenRx,
  MOSyncRx,
  MOSyncTx,
  numMetaOp = 3
};

static void
sendCollSync(struct fileFunnelStdio *bfd, int fileID, MPI_Comm commPio, int nProcsColl, MPI_Request (*req)[nProcsColl],
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
funnelFileCleanup(struct fileFunnelStdio *openFileFunnels, int fileID, size_t *openFileFunnelsFill, size_t *maxOpenFileIDp1,
                  size_t *numOpenRequests, MPI_Comm commPio, int nProcsColl)
{
  xdebug("all are finished with file %d, delete file table entry", fileID);
  struct fileFunnelStdio *bfd = openFileFunnels + fileID;
  assert(bfd->syncOffset == bfd->offset);
  int retval = destroyBFiledataP(bfd);
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

void
pioWriterStdIO(void)
{
  size_t openFileFunnelsSize = 1, openFileFunnelsFill = 0;
  struct fileFunnelStdio *openFileFunnels = Malloc(sizeof(*openFileFunnels) * openFileFunnelsSize);
  for (size_t i = 0; i < openFileFunnelsSize; ++i) openFileFunnels[i].name = NULL;
  MPI_Comm commPio = commInqCommPio();
  int nProcsColl = commInqSizeColl();

  int outstandingFinalizations = nProcsColl;
  struct cdiPioConf *conf = cdiPioGetConf();

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
  size_t maxOpenFileIDp1 = 0, numOpenSends = 0, numOpenRequests = (size_t) nProcsColl * (numMetaOp + maxOpenFileIDp1);
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
                struct fileFunnelStdio *bfd = openFileFunnels + fileID;
                if (!bfd->name)
                  {
                    for (size_t i = 0; i < prevOpenFileFunnelsSize; ++i)
                      if (openFileFunnels[i].name && !strcmp(openFileFunnels[i].name, filename))
                        xabort("filename %s is already open!"
                               " CDI-PIO does not support concurrent access"
                               " through different filehandles.",
                               filename);
                    initBFiledataP(bfd, filename, buffersize, (size_t) nProcsColl);
                    ++openFileFunnelsFill;
                  }
                else if (strcmp(filename, bfd->name) != 0)
                  xabort("filename is not consistent, fileID=%d,"
                         " \"%s\" vs. \"%s\"",
                         fileID, filename, bfd->name);
                int tag = encodeFileOpTag(fileID, IO_Send_buffer);
                xmpi(MPI_Irecv(bfd->db[source].buffer, (int) bfd->db[source].size, MPI_UNSIGNED_CHAR, source, tag, commPio,
                               req[numMetaOp + fileID] + source));
                xmpi(MPI_Irecv(reqBuffer, (int) openMsgSize, MPI_PACKED, source, IO_Open_file, commPio, req[0] + source));
                numOpenRequests = (size_t) nProcsColl * (numMetaOp + maxOpenFileIDp1);
              }
              break;

            case IO_Send_buffer:
              {
                int fileID = rcvIdx / nProcsColl - numMetaOp;
                assert(fileID >= 0 && (size_t) fileID < openFileFunnelsSize && openFileFunnels[fileID].name);
                struct fileFunnelStdio *bfd = openFileFunnels + fileID;

                xdebug("command %s, fileID=%d, name=%s", cdiPioCmdStrTab[command], fileID, bfd->name);

                int messagesize;
                xmpi(MPI_Get_count(statui + cmpltIdx, MPI_UNSIGNED_CHAR, &messagesize));
                size_t amount = (size_t) messagesize;

                writeP(bfd, source, amount);
                if (!bfd->outstandingSync && bfd->syncOffset == bfd->offset)
                  sendCollSync(bfd, fileID, commPio, nProcsColl, req, &numOpenSends, wstats);
                off_t *restrict perCollAmounts = bfd->perCollAmounts;
                /* end of collector stream not yet reached? */
                if (perCollAmounts[source] != perCollAmounts[nProcsColl + source])
                  {
                    int tag = encodeFileOpTag(fileID, IO_Send_buffer);
                    xmpi(MPI_Irecv(bfd->db[source].buffer, (int) bfd->db[source].size, MPI_UNSIGNED_CHAR, source, tag, commPio,
                                   req[numMetaOp + fileID] + source));
                  }
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
                struct fileFunnelStdio *bfd = openFileFunnels + fileID;
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
                struct fileFunnelStdio *bfd = openFileFunnels + fileID;
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

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
