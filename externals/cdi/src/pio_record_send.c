#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <unistd.h>

#include "cdipio.h"
#include "dmemory.h"
#include "namespace.h"
#include "pio.h"
#include "pio_comm.h"
#include "pio_dbuffer.h"
#include "pio_id_set.h"
#include "pio_impl.h"
#include "pio_util.h"

enum
{
  numRecordSendBuf = 2,
  numFlushOp = 2,
};
struct remoteFileBuf
{
  char *name;
  struct dBuffer dbuf[numRecordSendBuf];
  /* extra request for sync signal to and go-on signal from writer */
  MPI_Request request[numRecordSendBuf + numFlushOp];
  struct syncMsg syncMsg;
  off_t amount;
  unsigned dbufIdx;
  int tsID;
};

static struct remoteFileBuf *openRemoteFiles;
static unsigned openRemoteFilesSize, openRemoteFilesFill;

static void
initRemoteFileBuf(struct remoteFileBuf *restrict rfile, const char *filename, size_t bs)
{

  xdebug("filename=%s, buffersize=%zu, in", filename, bs);

  size_t len = strlen(filename);
  rfile->name = Malloc(len + 1);
  memcpy(rfile->name, filename, len + 1);
  rfile->tsID = 0;
  rfile->amount = 0;
  /* init output buffer */

  xdebug("filename=%s, init output buffer", filename);

  for (size_t i = 0; i < numRecordSendBuf; ++i) cdiPioDbufferInit(rfile->dbuf + i, bs);
  for (size_t i = 0; i < numRecordSendBuf + numFlushOp; ++i) rfile->request[i] = MPI_REQUEST_NULL;

  rfile->dbufIdx = 0;

  xdebug("added name=%s, return", rfile->name);
}

static int
destroyRemoteFileBuf(struct remoteFileBuf *restrict rfile)
{
  MPI_Status status[numRecordSendBuf + numFlushOp];

  xdebug("filename=%s, cleanup, in", rfile->name);

  xmpiStats(MPI_Waitall(numRecordSendBuf + numFlushOp, rfile->request, status), numRecordSendBuf, status);

  for (size_t i = 0; i < numRecordSendBuf; ++i) cdiPioDbufferDestroy(rfile->dbuf + i);
  Free(rfile->name);
  rfile->name = NULL;

  xdebug("%s", "cleaned up, return");

  return 0;
}

/***************************************************************/
/* send buffer to writer and swap buffer for filling */
static void
sendP(struct remoteFileBuf *rfile)
{
  MPI_Status statui[numFlushOp];
  if (rfile->request[numRecordSendBuf + 1] != MPI_REQUEST_NULL)
    xmpiStats(MPI_Waitall(numFlushOp, rfile->request + numRecordSendBuf, statui), numFlushOp, statui);
  size_t dbufIdx = rfile->dbufIdx;
  struct dBuffer *restrict dbuf = rfile->dbuf;
  size_t amount = cdiPioDbufferGetPos(dbuf + dbufIdx);
  if (amount)
    {
      int fileID = (int) (rfile - openRemoteFiles);
      int tag = encodeFileOpTag(fileID, IO_Send_buffer);

      xdebug("send buffer for %s, size: %zu bytes, command=%s, in", rfile->name, amount, cdiPioCmdStrTab[IO_Send_buffer]);

      /* FIXME: amount > INT_MAX unhandled */
      xmpi(MPI_Isend(dbuf[dbufIdx].buffer, (int) amount, MPI_UNSIGNED_CHAR, commInqSizePio() - 1, tag, commInqCommPio(),
                     rfile->request + dbufIdx));
      rfile->amount += (off_t) amount;

      /* change outputBuffer */
      dbufIdx ^= 1;
      xmpiStat(MPI_Wait(rfile->request + dbufIdx, statui), statui);
      cdiPioDbufferReset(dbuf + dbufIdx);
      rfile->dbufIdx = (unsigned) dbufIdx;
    }
}

static void
flushOp(struct remoteFileBuf *rfile, int tsID)
{
  assert(rfile != NULL && tsID == rfile->tsID + 1);
  sendP(rfile);
  rfile->tsID = tsID;
  int fileID = (int) (rfile - openRemoteFiles), funnelRank = commInqSizePio() - 1;
  MPI_Comm commPio = commInqCommPio();
  MPI_Request *req = rfile->request + numRecordSendBuf;
  rfile->syncMsg = (struct syncMsg){ .amount = rfile->amount, .fileID = fileID, .command = IO_Sync_file };
  xmpi(MPI_Isend(&rfile->syncMsg, 1, cdiPioSyncMsgDt, funnelRank, IO_Sync_file, commPio, req));
  int syncFromFunnelTag = encodeFileOpTag(fileID, IO_Sync_file);
  xmpi(MPI_Irecv(NULL, 0, MPI_INT, funnelRank, syncFromFunnelTag, commPio, req + 1));
  MPI_Status statui[numFlushOp];
  int numComplete, completeIdx[numFlushOp];
  xmpiStats(MPI_Testsome(numFlushOp, req, &numComplete, completeIdx, statui), numComplete, statui);
}

static size_t
pioSendWrite(int fileID, const void *buffer, size_t len, int tsID)
{
  assert(fileID >= 0 && (size_t) fileID < openRemoteFilesSize && openRemoteFiles[fileID].name);
  struct remoteFileBuf *rfile = openRemoteFiles + fileID;

  bool flush = tsID != rfile->tsID;

  size_t dbufIdx = rfile->dbufIdx;
  struct dBuffer *restrict dbuf = rfile->dbuf;

  if (flush)
    {
      xdebug("tsID = %d, flush buffer for fileID=%d", tsID, fileID);
      flushOp(rfile, tsID);
      dbufIdx ^= 1;
    }

  int filled = cdiPioDbufferAppend(dbuf + dbufIdx, buffer, len);

  xdebug("id = %d, tsID = %d, pushed %lu byte data on buffer, filled = %d", fileID, tsID, len, filled);

  int error = 0;
  if (filled == 1)
    {
      if (flush)
        error = filled;
      else
        {
          sendP(rfile);
          dbufIdx = rfile->dbufIdx;
          error = cdiPioDbufferAppend(dbuf + dbufIdx, buffer, len);
        }
    }

  if (error == 1) xabort("did not succeed filling output buffer, id=%d", fileID);

  return len;
}

static int
pioSendClose(int fileID)
{
  assert(fileID >= 0 && (size_t) fileID < openRemoteFilesSize && openRemoteFiles[fileID].name);

  struct remoteFileBuf *rfile = openRemoteFiles + fileID;

  xdebug("fileID %d: send buffer, close file and cleanup", fileID);

  sendP(rfile);
  MPI_Comm commPio = commInqCommPio();
  int funnelRank = commInqSizePio() - 1;
  struct syncMsg query = { .amount = rfile->amount, .fileID = fileID, .command = IO_Close_file };
  xmpi(MPI_Send(&query, 1, cdiPioSyncMsgDt, funnelRank, IO_Sync_file, commPio));
  /* remove file element and wait for pending data sends to finish */
  destroyRemoteFileBuf(rfile);
  /* wait for other collectors to also close the file
   * this prevents one process from re-using the file ID before
   * another has sent the close */
  int iret;
  xmpi(MPI_Bcast(&iret, 1, MPI_INT, funnelRank, commPio));

  return iret;
}

static int
pioSendOpen(const char *filename, const char *mode)
{
  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  struct cdiPioConf *conf = cdiPioGetConf();

  /* init and add struct remoteFileBuf */
  for (size_t i = 0; i < openRemoteFilesSize; ++i)
    if (openRemoteFiles[i].name && !strcmp(openRemoteFiles[i].name, filename))
      {
        Warning("filename %s is already open!"
                " CDI-PIO does not support concurrent access"
                " through different filehandles.",
                filename);
        return CDI_EINVAL;
      }

  int fileID = CDI_ELIMIT;
  if (openRemoteFilesSize == openRemoteFilesFill)
    {
      fileID = (int) openRemoteFilesSize;
      if (openRemoteFilesSize == (size_t) INT_MAX + 1) return CDI_ELIMIT;
      openRemoteFilesSize = openRemoteFilesSize ? openRemoteFilesSize * 2 : 4;
      if (openRemoteFilesSize > (size_t) INT_MAX + 1) openRemoteFilesSize = (size_t) INT_MAX + 1;
      openRemoteFiles = Realloc(openRemoteFiles, sizeof(*openRemoteFiles) * openRemoteFilesSize);
      for (size_t i = (size_t) fileID; i < openRemoteFilesSize; ++i) openRemoteFiles[i].name = NULL;
    }
  else
    {
      for (size_t i = 0; i < openRemoteFilesSize; ++i)
        if (openRemoteFiles[i].name == NULL)
          {
            fileID = (int) i;
            break;
          }
    }
  struct remoteFileBuf *rfile = openRemoteFiles + fileID;
  initRemoteFileBuf(rfile, filename, conf->writeAggBufLim);

  xdebug("filename=%s, init and added remoteFileBuf, return fileID = %d", filename, fileID);

  /* put filename, id and buffersize on buffer */
  size_t dbufIdx = rfile->dbufIdx;
  MPI_Comm commPio = commInqCommPio();
  int intPackSize;
  xmpi(MPI_Pack_size(2, MPI_INT, commPio, &intPackSize));
  size_t maxPathLen = (size_t) conf->maxPathLen, openMsgSize = (size_t) intPackSize + maxPathLen;

  unsigned char *buf = rfile->dbuf[dbufIdx].size >= openMsgSize ? rfile->dbuf[dbufIdx].buffer : Malloc(openMsgSize);
  size_t bufSize = conf->writeAggBufLim;
  size_t nameLen = strlen(filename);
  assert(nameLen <= (size_t) conf->maxPathLen);
  int position = 0;
  int openMsgHdr[2] = { [0] = fileID, [1] = (int) nameLen };
  int funnelRank = commInqSizePio() - 1;
  xmpi(MPI_Pack(openMsgHdr, 2, MPI_INT, buf, (int) bufSize, &position, commPio));
  xmpi(MPI_Pack((void *) filename, (int) nameLen, MPI_CHAR, buf, (int) bufSize, &position, commPio));
  xmpi(MPI_Send(buf, position, MPI_PACKED, funnelRank, IO_Open_file, commPio));
  if (buf != rfile->dbuf[dbufIdx].buffer) Free(buf);
  return fileID;
}

static void
pioSendFinalize(void)
{
  int funnelRank = commInqSizePio() - 1;
  MPI_Comm commPio = commInqCommPio();

  static const struct syncMsg finQuery = { .amount = 0, .fileID = 0, .command = IO_Finalize };
  xmpi(MPI_Send((void *) &finQuery, 1, cdiPioSyncMsgDt, funnelRank, IO_Sync_file, commPio));

  xdebug("%s", "SENT MESSAGE WITH TAG \"IO_FINALIZE\" TO SPECIAL PROCESS");

  if (openRemoteFilesFill)
    xabort("files still open.");
  else
    {
      xdebug("%s", "destroy set");
      Free(openRemoteFiles);
    }
  cdiPioDestroySyncMsgDt();
}

void
pioSendInitialize(void)
{
  int numIOServers = commInqSizePio();
  cdiPioCreateSyncMsgDt();
  if (numIOServers < 2) xabort("error: # of I/O processes must be >= 2 for mode, but is %d", numIOServers);

  int isCollector = commInqRankColl() != -1;
  if (!isCollector)
    {
      switch (commInqIOMode())
        {
        case PIO_WRITER: pioWriterStdIO(); break;
#ifdef _POSIX_ASYNCHRONOUS_IO
        case PIO_ASYNCH: pioWriterAIO(); break;
#endif
        }
      cdiPioDestroySyncMsgDt();
    }
  else
    {
      namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(pioSendOpen));
      namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(pioSendClose));
      namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(pioSendWrite));
      namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(pioSendFinalize));
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
