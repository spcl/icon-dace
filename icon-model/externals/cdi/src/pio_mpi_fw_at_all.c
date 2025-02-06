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
#include "pio_impl.h"
#include "pio_util.h"

struct fileMPIFWAA
{
  char *name;
  MPI_File fh;
  int fileID;
  MPI_Offset pos;
  int *collWriteSize;
};

static struct fileMPIFWAA *openFiles;
static unsigned openFilesSize, openFilesFill;

/***************************************************************/

static void
initAFiledataFileWriteAtAll(struct fileMPIFWAA *of, const char *filename, size_t bs)
{
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio();
  size_t nameSize = strlen(filename) + 1;
  of->collWriteSize = Malloc(sizeof(of->collWriteSize[0]) * (size_t) sizePio);
  of->name = Malloc(nameSize);
  memcpy(of->name, filename, nameSize);

  MPI_Info open_info = MPI_INFO_NULL;
  xmpi(MPI_Info_create(&open_info));
  xmpi(MPI_Info_set(open_info, "access_style", "write_once"));
  xmpi(MPI_Info_set(open_info, "collective_buffering", "true"));
  /* tell IBM PE to buffer just as much as one buffer holds */
  {
    char buf_size_str[3 * sizeof(size_t) * CHAR_BIT / 8 + 1];
    snprintf(buf_size_str, sizeof(buf_size_str), "%zu", bs);
    xmpi(MPI_Info_set(open_info, "IBM_io_buffer_size", buf_size_str));
    xmpi(MPI_Info_set(open_info, "IBM_largeblock_io", "false"));
  }
  xmpi(MPI_File_open(commPio, of->name, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN, open_info, &of->fh));
  xmpi(MPI_Info_free(&open_info));
  of->pos = 0;
}

/***************************************************************/

static int
destroyAFiledataFileWriteAtAll(struct fileMPIFWAA *of)
{
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
  Free(of->collWriteSize);
  Free(of->name);
  of->name = NULL;

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static size_t
fwFileWriteAtAll(int fileID, const void *buffer, size_t len)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWAA *of = openFiles + fileID;
  xassert(len <= INT_MAX);
  MPI_Comm commPio = commInqCommPio();
  int sizePio = commInqSizePio(), rankPio = commInqRankPio();
  /* find position to write to */
  of->collWriteSize[rankPio] = (int) len;
  xmpi(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, of->collWriteSize, 1, MPI_INT, commPio));
  MPI_Offset myPos = of->pos, nextWritePos;
  for (size_t i = 0; i < (size_t) rankPio; ++i) myPos += of->collWriteSize[i];
  nextWritePos = myPos;
  for (size_t i = (size_t) rankPio; i < (size_t) sizePio; ++i) nextWritePos += of->collWriteSize[i];
  /* write buffer */
  xmpi(MPI_File_write_at_all(of->fh, myPos, (void *) buffer, (int) len, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE));
  of->pos = nextWritePos;
  return len;
}

/***************************************************************/

static int
fcFileWriteAtAll(int fileID)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWAA *of = openFiles + fileID;
  int iret = destroyAFiledataFileWriteAtAll(of);
  --openFilesFill;
  return iret;
}

/***************************************************************/

static int
fowFileWriteAtAll(const char *filename, const char *mode)
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
  struct fileMPIFWAA *of = openFiles + fileID;
  ++openFilesFill;
  initAFiledataFileWriteAtAll(of, filename, conf->writeAggBufLim);
  return (int) fileID;
}

/***************************************************************/

static void
finalizeFileWriteAtAll(void)
{
  if (openFilesFill)
    xabort("files still open on exit!");
  else
    Free(openFiles);
}

/***************************************************************/

void
cdiPioFileWriteAtAllInit(void)
{
  namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(fowFileWriteAtAll));
  namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(fcFileWriteAtAll));
  namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(fwFileWriteAtAll));
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(finalizeFileWriteAtAll));
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
