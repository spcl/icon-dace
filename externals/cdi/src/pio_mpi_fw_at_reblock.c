#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <mpi.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "error.h"
#include "namespace.h"
#include "pio.h"
#include "pio_comm.h"
#include "pio_conf.h"
#include "pio_impl.h"
#include "pio_rpc.h"
#include "pio_util.h"

enum direction
{
  INACTIVE_TX,
  RECV,
  SEND,
};

struct IOmsg
{
  int pos, len; /* local index of block buffer
                 * affected (RECV only) */
  void *sendBuf;
  int rank;
  enum direction direction;
};

struct pendingBufWrite
{
  int pass;
  int incoming;
};

struct fileMPIFWAR
{
  MPI_File fh;
  int numMsg;
  MPI_Offset pos;
  char *name;
  int blockSize;
  int numBlockBuf;                 /* number of block-sized buffers
                                    * pointed to by blockBuf */
  struct pendingBufWrite *pending; /* pendingBufWrite[i].pass
                                    * indicates pass if non-negative,
                                    * count gives number of
                                    * outstanding recvs still to be
                                    * processed for block buf i */
  unsigned char *blockBuf;
  int msgSize;
  struct IOmsg *msgs;
  MPI_Request *reqs;
  long *collWriteSize; /* used to allgather sizes of writes
                        * on different processes */
};

static struct fileMPIFWAR *openFiles;
static unsigned openFilesSize, openFilesFill;

/***************************************************************/

static inline void
initReblockPendingMsg(struct fileMPIFWAR *of, size_t i)
{
  of->msgs[i].pos = -1;
  of->msgs[i].len = -1;
  of->msgs[i].sendBuf = NULL;
  of->msgs[i].rank = -1;
  of->msgs[i].direction = INACTIVE_TX;
  of->reqs[i] = MPI_REQUEST_NULL;
}

static inline size_t
gcd(size_t a, size_t b)
{
  while (b != 0)
    {
      size_t t = a;
      a = b;
      b = t % b;
    }
  return a;
}

static inline size_t
lcm(size_t a, size_t b)
{
  return a / gcd(a, b) * b;
}

static size_t
getXferBufAlign(const char *path)
{
  long align = -1L;
#if HAVE_DECL__PC_REC_XFER_ALIGN
  align = pathconf(path, _PC_REC_XFER_ALIGN);
#endif
  if (align == -1L)
    align =
#if HAVE_DECL_POSIX_REC_XFER_ALIGN
        POSIX_REC_XFER_ALIGN
#else
        commonPageSize
#endif
        ;
  return (size_t) align;
}

static void
initAFiledataFileWriteAtReblock(struct fileMPIFWAR *of, const char *filename, const struct cdiPioConf *conf)
{
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio();
  size_t nameSize = strlen(filename) + 1;

  MPI_Info open_info = MPI_INFO_NULL;
  xmpi(MPI_Info_create(&open_info));
  xmpi(MPI_Info_set(open_info, "access_style", "write_once"));
  xmpi(MPI_Info_set(open_info, "collective_buffering", "false"));
  /* tell IBM PE to not buffer anything, we block-align all writes */
  xmpi(MPI_Info_set(open_info, "IBM_io_buffer_size", "0"));
  xmpi(MPI_Info_set(open_info, "IBM_largeblock_io", "true"));
  MPI_File fh;
  xmpi(MPI_File_open(commPio, (char *) filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN, open_info, &fh));
  xmpi(MPI_Info_free(&open_info));
  /* find block size of underlying file system */
  size_t blockSize, bufBlockAlign;
  {
    struct stat posixstat;
    int rc = stat(filename, &posixstat);
    if (rc < 0)
      {
        perror("failed to stat file after open, block size unavailable,"
               " assuming 4MiB");
        posixstat.st_blksize = 4 * 1024 * 1024;
      }
    blockSize = (size_t) posixstat.st_blksize < SIZE_MAX ? (size_t) posixstat.st_blksize : SIZE_MAX / 2;
    /* prevent inefficiently small writes on granular file systems */
#define MINBLOCKSIZE ((size_t) 1 << 19)
    if ((blockSize < MINBLOCKSIZE) & !(MINBLOCKSIZE % blockSize)) blockSize = MINBLOCKSIZE;
#undef MINBLOCKSIZE
    /* ensure block size also meets page and direct I/O buffer
     * alignment requirements */
    size_t pageSize = cdiGetPageSize(conf->largePageAlign), IOAlign = getXferBufAlign(filename);
    bufBlockAlign = lcm(pageSize, IOAlign);
    blockSize = lcm(blockSize, bufBlockAlign);
  }
  /* round bufSize to next multiple of block size and I/O alignment */
  size_t bufSize = roundUpToMultiple(conf->writeAggBufLim, blockSize);
  size_t numBlockBuf = bufSize / blockSize;
  assert(blockSize <= INT_MAX && numBlockBuf <= INT_MAX);
  /* never go less than double-buffering */
  if (numBlockBuf < 2)
    {
      numBlockBuf = 2;
      bufSize = blockSize * 2;
    }
  of->collWriteSize = Malloc(sizeof(of->collWriteSize[0]) * (size_t) sizePio);
  of->fh = fh;
  of->name = Malloc(nameSize);
  memcpy(of->name, filename, nameSize);
  {
    void *ptr = NULL;
    int err = posix_memalign(&ptr, bufBlockAlign, bufSize);
    if (!err)
      of->blockBuf = ptr;
    else
      cdiAbort(__FILE__, __func__, __LINE__, "posix_memalign failed: %s", strerror(errno));
  }
  of->pending = (struct pendingBufWrite *) Malloc(numBlockBuf * sizeof(struct pendingBufWrite));
  for (size_t i = 0; i < numBlockBuf; ++i)
    {
      of->pending[i].pass = -1;
      of->pending[i].incoming = 0;
    }
  of->blockSize = (int) blockSize;
  of->numBlockBuf = (int) numBlockBuf;
  size_t numPeers = (size_t) commInqSizePio();
  /* start with 2 sends and 2 recvs per peer simultaneously */
  size_t msgSize = numPeers * 4;
  assert(msgSize <= INT_MAX);
  of->msgSize = (int) msgSize;
  of->msgs = Malloc(sizeof(of->msgs[0]) * msgSize);
  of->reqs = Malloc(sizeof(of->reqs[0]) * msgSize);
  for (size_t i = 0; i < (size_t) msgSize; ++i) initReblockPendingMsg(of, i);
  of->pos = 0;
  of->numMsg = 0;
}

/***************************************************************/

static void flushReblockBuffer(struct fileMPIFWAR *of, int blockBufIdx);

static int
destroyAFiledataFileWriteAtReblock(struct fileMPIFWAR *restrict of)
{
  size_t numBlockBuf = (size_t) of->numBlockBuf;
  /* flush pending buffers */
  /** 1. handle all outstanding messages */
  xmpi(MPI_Waitall(of->numMsg, of->reqs, MPI_STATUSES_IGNORE));
  of->numMsg = 0;
  for (size_t block = 0; block < numBlockBuf; ++block)
    {
      of->pending[block].incoming = 0;
      if (of->pending[block].pass != -1) flushReblockBuffer(of, (int) block);
    }
  /* close file */
  MPI_Offset endpos, fsize;
  endpos = of->pos;
  xmpi(MPI_File_get_size(of->fh, &fsize));
  /* does the file need to be truncated? */
  MPI_Comm commPio = commInqCommPio();
  int trailingOctets = fsize > endpos;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &trailingOctets, 1, MPI_INT, MPI_LOR, commPio));
  if (trailingOctets) xmpi(MPI_File_set_size(of->fh, endpos));
  int iret = MPI_File_close(&of->fh);
  for (size_t i = 0; i < (size_t) of->msgSize; ++i) Free(of->msgs[i].sendBuf);
  free(of->blockBuf);
  Free(of->msgs);
  Free(of->reqs);
  Free(of->pending);
  Free(of->collWriteSize);
  Free(of->name);
  of->name = NULL;

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static inline long
lmin(long a, long b)
{
  return a < b ? a : b;
}

static void
flushReblockBuffer(struct fileMPIFWAR *of, int blockBufIdx)
{
  int blockSize = of->blockSize;
  unsigned char *blockBuf = of->blockBuf + blockSize * blockBufIdx;
  int finIdx, numMsg = of->numMsg;
  while (of->pending[blockBufIdx].incoming)
    {
      /* todo: switch to MPI_Waitsome */
      xmpi(MPI_Waitany(numMsg, of->reqs, &finIdx, MPI_STATUS_IGNORE));
      if (finIdx != MPI_UNDEFINED)
        {
          int blockBufIdxOfTx = of->msgs[finIdx].pos;
          if (blockBufIdxOfTx >= 0 && of->msgs[finIdx].direction == RECV)
            --(of->pending[blockBufIdxOfTx].incoming);
          else if (of->msgs[finIdx].direction == SEND)
            Free(of->msgs[finIdx].sendBuf);
          else
            xabort("internal error");
          of->msgs[finIdx] = of->msgs[numMsg - 1];
          of->msgs[numMsg - 1].sendBuf = NULL;
          of->reqs[finIdx] = of->reqs[numMsg - 1];
          --numMsg;
        }
    }
  int sizePio = commInqSizePio(), rankPio = commInqRankPio();
  MPI_Offset ofs = (MPI_Offset) blockSize
                   * (((MPI_Offset) of->pending[blockBufIdx].pass * (MPI_Offset) sizePio * of->numBlockBuf)
                      + (MPI_Offset) blockBufIdx * (MPI_Offset) sizePio + (MPI_Offset) rankPio);
  int wsize = (of->pos >= ofs + blockSize) ? blockSize : (int) (of->pos - ofs);
  xmpi(MPI_File_write_at(of->fh, ofs, blockBuf, wsize, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE));
  of->pending[blockBufIdx].pass = -1;
  of->numMsg = numMsg;
}

static void
reblockMoreMsgs(struct fileMPIFWAR *of, int numMsg)
{
  /* optimize with MPI_Testsome */
  if (of->msgSize == numMsg)
    {
      size_t newMsgSize = (size_t) numMsg * 2;
      of->msgs = Realloc(of->msgs, sizeof(of->msgs[0]) * newMsgSize);
      of->reqs = Realloc(of->reqs, sizeof(of->reqs[0]) * newMsgSize);
      for (size_t i = (size_t) numMsg; i < newMsgSize; ++i) initReblockPendingMsg(of, i);
      assert(newMsgSize <= INT_MAX);
      of->msgSize = (int) newMsgSize;
    }
}

static size_t
fwFileWriteAtReblock(int fileID, const void *buffer, size_t len)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWAR *of = openFiles + fileID;
  xassert(len <= INT_MAX);
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio(), rankPio = commInqRankPio();
  /* find position to write to */
  of->collWriteSize[rankPio] = (long) len;
  xmpi(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, of->collWriteSize, 1, MPI_LONG, commPio));
  /* figure out which block buffers intersect locally held data and
   * what remotely held data intersects buffers on task */
  int blockSize = of->blockSize;
  MPI_Offset fWOfs = of->pos;
  int numBlockBuf = of->numBlockBuf;
  int numMsg = of->numMsg;
  const unsigned char *inBuf = buffer;
  const unsigned char *directWriteBuf = NULL;
  /* positive value iff direct write occurs */
  int directWriteSize = -1;
  MPI_Offset directWriteOfs = 0;
  for (int collRank = 0; collRank < sizePio; ++collRank)
    if (of->collWriteSize[collRank])
      {
        long remaining = of->collWriteSize[collRank];
        do
          {
            MPI_Offset collBlockSize = (MPI_Offset) blockSize * (MPI_Offset) sizePio;
            int pass = (int) (fWOfs / collBlockSize / numBlockBuf);
            int inBlockPos = (int) (fWOfs % (MPI_Offset) blockSize);
            int txLen = (int) lmin(blockSize - inBlockPos, remaining);
            int destRank = (int) (fWOfs / (MPI_Offset) blockSize) % sizePio;
            if (txLen == blockSize)
              {
                /* a properly aligned portion can be written directly */
                txLen = (int) (remaining - remaining % blockSize);
                if (collRank == rankPio)
                  {
                    directWriteBuf = inBuf;
                    inBuf += txLen;
                    directWriteSize = txLen;
                    directWriteOfs = fWOfs;
                  }
              }
            else if (destRank == rankPio)
              {
                int blockBufIdx = (int) (fWOfs / collBlockSize) % numBlockBuf;
                if (of->pending[blockBufIdx].pass >= 0 && of->pending[blockBufIdx].pass != pass)
                  {
                    of->numMsg = numMsg;
                    flushReblockBuffer(of, blockBufIdx);
                    numMsg = of->numMsg;
                  }
                if (collRank != rankPio)
                  {
                    reblockMoreMsgs(of, numMsg);
                    /* this rank is to write out (part of) the data, but it
                     * resides on another rank */
                    of->msgs[numMsg] = (struct IOmsg){
                      .pos = blockBufIdx,
                      .len = txLen,
                      .rank = collRank,
                      .direction = RECV,
                    };
                    xmpi(MPI_Irecv(of->blockBuf + (size_t) blockBufIdx * (size_t) blockSize + (size_t) inBlockPos, txLen,
                                   MPI_UNSIGNED_CHAR, collRank, BLOCK_XFER, commPio, of->reqs + numMsg));
                    ++numMsg;
                    ++(of->pending[blockBufIdx].incoming);
                  }
                else /* if (collRank == rankPio) */
                  {
                    memcpy(of->blockBuf + (size_t) blockBufIdx * (size_t) blockSize + (size_t) inBlockPos, inBuf, (size_t) txLen);
                    inBuf += txLen;
                  }
                of->pending[blockBufIdx].pass = pass;
              }
            else if (collRank == rankPio)
              {
                reblockMoreMsgs(of, numMsg);
                /* this rank has the data and will send it to the one doing
                 * the writing */
                of->msgs[numMsg] = (struct IOmsg){
                  .pos = -1,
                  .len = txLen,
                  .rank = destRank,
                  .direction = SEND,
                };
                void *restrict buf = of->msgs[numMsg].sendBuf = Realloc(of->msgs[numMsg].sendBuf, (size_t) txLen);
                memcpy(buf, inBuf, (size_t) txLen);
                xmpi(MPI_Isend(buf, txLen, MPI_UNSIGNED_CHAR, destRank, BLOCK_XFER, commPio, of->reqs + numMsg));
                inBuf += txLen;
                ++numMsg;
              }
            fWOfs += txLen;
            remaining -= txLen;
          }
        while (remaining);
      }
  if (directWriteSize > -1)
    xmpi(MPI_File_write_at(of->fh, directWriteOfs, (unsigned char *) directWriteBuf, directWriteSize, MPI_UNSIGNED_CHAR,
                           MPI_STATUS_IGNORE));
  of->numMsg = numMsg;
  of->pos = fWOfs;
  return len;
}

/***************************************************************/

static int
fcFileWriteAtReblock(int fileID)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWAR *of = openFiles + fileID;
  int iret = destroyAFiledataFileWriteAtReblock(of);
  --openFilesFill;
  return iret;
}

/***************************************************************/

static int
fowFileWriteAtReblock(const char *filename, const char *mode)
{
  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  const struct cdiPioConf *conf = cdiPioGetConf();

  for (size_t i = 0; i < openFilesSize; ++i)
    if (openFiles[i].name && !strcmp(openFiles[i].name, filename))
      {
        Warning("filename %s is already open!"
                " CDI-PIO does not support concurrent access"
                " through different filehandles.",
                filename);
        return CDI_EINVAL;
      }

  size_t fileID = SIZE_MAX;
  if (openFilesSize == openFilesFill)
    {
      fileID = openFilesSize;
      if (openFilesSize == (size_t) INT_MAX + 1) return CDI_ELIMIT;
      openFilesSize = openFilesSize ? openFilesSize * 2 : 4;
      if (openFilesSize > (size_t) INT_MAX + 1) openFilesSize = (size_t) INT_MAX + 1;
      openFiles = Realloc(openFiles, sizeof(*openFiles) * openFilesSize);
      for (size_t i = fileID; i < openFilesSize; ++i) openFiles[i].name = NULL;
    }
  else
    {
      for (size_t i = 0; i < openFilesSize; ++i)
        if (openFiles[i].name == NULL)
          {
            fileID = i;
            break;
          }
    }
  struct fileMPIFWAR *of = openFiles + fileID;
  ++openFilesFill;
  initAFiledataFileWriteAtReblock(of, filename, conf);
  return (int) fileID;
}

/***************************************************************/

static void
finalizeFileWriteAtReblock(void)
{
  if (openFilesFill)
    xabort("files still open on exit!");
  else
    Free(openFiles);
}

/***************************************************************/

void
cdiPioFileWriteAtReblockInit(void)
{
  namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(fowFileWriteAtReblock));
  namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(fcFileWriteAtReblock));
  namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(fwFileWriteAtReblock));
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(finalizeFileWriteAtReblock));
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
