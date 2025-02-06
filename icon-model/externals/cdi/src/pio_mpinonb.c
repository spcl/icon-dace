#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "cdi.h"
#include "dmemory.h"
#include "namespace.h"
#include "pio.h"
#include "pio_comm.h"
#include "pio_dbuffer.h"
#include "pio_impl.h"
#include "pio_util.h"

struct fileMPIFWS
{
  char *name;
  struct dBuffer db[2];
  MPI_File fh;
  MPI_Request request;
  int tsID;
  unsigned dbufIdx;
};

static struct fileMPIFWS *openFiles;
static unsigned openFilesSize, openFilesFill;

/* Open MPI 2.0.2 forgets to set request to MPI_REQUEST_NULL */
#ifdef OMPI_MINOR_VERSION
#if OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION == 2
#define POSTWAIT_WORKAROUND(req) req = MPI_REQUEST_NULL
#else
#define POSTWAIT_WORKAROUND(req)
#endif
#else
#define POSTWAIT_WORKAROUND(req)
#endif

/***************************************************************/

static void
initAFiledataMPINONB(struct fileMPIFWS *of, const char *filename, size_t bufSize)
{
  MPI_Comm commPio = commInqCommPio();
  {
    size_t nameSize = strlen(filename) + 1;
    char *name = of->name = Malloc(nameSize);
    memcpy(name, filename, nameSize);
  }

  /* init output buffer */
  for (size_t dbufIdx = 0; dbufIdx < 2; ++dbufIdx) cdiPioDbufferInit(of->db + dbufIdx, bufSize);

  of->dbufIdx = 0;

  of->tsID = -1;

  MPI_Info open_info = MPI_INFO_NULL;
  /* tell IBM PE to buffer just as much as one buffer holds */
  {
    xmpi(MPI_Info_create(&open_info));
    char buf_size_str[3 * sizeof(size_t) * CHAR_BIT / 8 + 1];
    snprintf(buf_size_str, sizeof(buf_size_str), "%zu", bufSize);
    xmpi(MPI_Info_set(open_info, "IBM_io_buffer_size", buf_size_str));
    xmpi(MPI_Info_set(open_info, "IBM_largeblock_io", "true"));
  }
  xmpi(MPI_File_open(commPio, of->name, MPI_MODE_CREATE | MPI_MODE_WRONLY, open_info, &of->fh));
  xmpi(MPI_Info_free(&open_info));
#ifdef OMPI_MINOR_VERSION
  /* Bull X MPI 1.2.x has a defective shared file pointer that can be
   * righted with a collective operation, thus we write zero bytes
   * from all tasks iff we detect that this is Bull X MPI, which
   * fortunately can be found from a mismatch in MPI and OMPI versions. */
#if OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION == 2 && MPI_SUBVERSION == 1
#warning "implementing Bull X MPI work-around"
  {
    MPI_Status status;
    xmpiStat(MPI_File_write_ordered(of->fh, filename, 0, MPI_BYTE, &status), &status);
  }
#endif
#endif
  of->request = MPI_REQUEST_NULL;
}

/***************************************************************/

static int
destroyAFiledataMPINONB(struct fileMPIFWS *of)
{
  int iret = 0;
  MPI_Status status;
  MPI_Offset endpos;

  xdebug("IOPE%d: close file %d, name=\"%s\"", commInqRankGlob(), (int) (of - openFiles), of->name);

  /* close file */
  xmpiStat(MPI_Wait(&of->request, &status), &status);
  xmpi(MPI_Barrier(commInqCommPio()));
  xmpi(MPI_File_get_position_shared(of->fh, &endpos));
  xmpi(MPI_File_set_size(of->fh, endpos));
  iret = MPI_File_close(&(of->fh));

  /* file closed, cleanup */

  cdiPioDbufferDestroy(of->db);
  cdiPioDbufferDestroy(of->db + 1);
  Free(of->name);
  of->name = NULL;

  xdebug("IOPE%d: closed file, cleaned up, return", commInqRankGlob());

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static void
writeMPINONB(struct fileMPIFWS *of)
{
  MPI_Status status;
  int fileID = (int) (of - openFiles);

  /* write buffer */
  size_t dbufIdx = of->dbufIdx;
  int amount = (int) cdiPioDbufferGetPos(of->db + dbufIdx);

  if (amount == 0) return;

  xdebug3("IOPI%d: Write buffer, size %d bytes, in", commInqRankGlob(), amount);

  xmpiStat(MPI_Wait(&of->request, &status), &status);
  xmpi(MPI_File_iwrite_shared(of->fh, of->db[dbufIdx].buffer, amount, MPI_UNSIGNED_CHAR, &of->request));
  xdebug("%d bytes written for fileID=%d", amount, fileID);

  /* change outputBuffer */

  cdiPioDbufferReset(of->db + dbufIdx);
  of->dbufIdx = (unsigned) (dbufIdx ^ 1);
}

/***************************************************************/

static size_t
fwMPINONB(int fileID, const void *buffer, size_t len, int tsID)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWS *of = openFiles + fileID;

  bool flush = tsID != of->tsID;

  if (flush)
    {
      xdebug3("IOPE%d: tsID = %d, flush buffer", commInqRankPio(), tsID);
      writeMPINONB(of);
      of->tsID = tsID;
      MPI_Status status;
      xmpiStat(MPI_Wait(&of->request, &status), &status);
      POSTWAIT_WORKAROUND(of->request);
      xmpi(MPI_Barrier(commInqCommPio()));
    }

  size_t dbufIdx = of->dbufIdx;
  int filled = cdiPioDbufferAppend(of->db + dbufIdx, buffer, len);

  xdebug3("IOPE%d: fileID = %d, tsID = %d,"
          " pushed data on buffer, filled = %d",
          commInqRankPio(), fileID, tsID, filled);

  int error = 0;
  if (filled == 1)
    {
      if (flush)
        error = filled;
      else
        {
          writeMPINONB(of);
          dbufIdx ^= 1;
          error = cdiPioDbufferAppend(of->db + dbufIdx, buffer, len);
        }
    }

  if (error == 1) xabort("did not succeed filling output buffer, fileID=%d", fileID);

  return len;
}

/***************************************************************/

static int
fcMPINONB(int fileID)
{
  xdebug("IOPE%d: write buffer, close file and cleanup, in %d", commInqRankPio(), fileID);

  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWS *of = openFiles + fileID;

  writeMPINONB(of);
  MPI_Status status;
  xmpiStat(MPI_Wait(&(of->request), &status), &status);
  POSTWAIT_WORKAROUND(of->request);
  /* remove file element */
  int iret = destroyAFiledataMPINONB(of);
  --openFilesFill;
  return iret;
}

/***************************************************************/

static int
fowMPINONB(const char *filename, const char *mode)
{
  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  struct cdiPioConf *conf = cdiPioGetConf();

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
  struct fileMPIFWS *of = openFiles + fileID;
  ++openFilesFill;
  initAFiledataMPINONB(of, filename, conf->writeAggBufLim);
  return (int) fileID;
}

/***************************************************************/

static void
finalizeMPINONB(void)
{
  if (openFilesFill) xabort("files still open on exit!");
  Free(openFiles);
}

/***************************************************************/

void
initMPINONB(void)
{
  namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(fowMPINONB));
  namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(fcMPINONB));
  namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(fwMPINONB));
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(finalizeMPINONB));
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
