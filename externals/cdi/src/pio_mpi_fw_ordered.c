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

struct fileMPIFWO
{
  char *name;
  MPI_File fh;
  int fileID;
};

static struct fileMPIFWO *openFiles;
static unsigned openFilesSize, openFilesFill;

/***************************************************************/

static void
initAFiledataFileWriteOrdered(struct fileMPIFWO *of, const char *filename, size_t bs)
{
  MPI_Comm commPio = commInqCommPio();
  size_t nameSize = strlen(filename) + 1;
  of->name = Malloc(nameSize);
  strcpy(of->name, filename);

  MPI_Info open_info = MPI_INFO_NULL;
  xmpi(MPI_Info_create(&open_info));
  xmpi(MPI_Info_set(open_info, "access_style", "sequential,write_once"));
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
}

/***************************************************************/

static int
destroyAFiledataFileWriteOrdered(struct fileMPIFWO *of)
{
  /* close file */
  MPI_Offset endpos, fsize;
  xmpi(MPI_File_get_position_shared(of->fh, &endpos));
  xmpi(MPI_File_get_size(of->fh, &fsize));
  /* does the file need to be truncated? */
  MPI_Comm commPio = commInqCommPio();
  int trailingOctets = fsize > endpos;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &trailingOctets, 1, MPI_INT, MPI_LOR, commPio));
  if (trailingOctets) xmpi(MPI_File_set_size(of->fh, endpos));
  int iret = MPI_File_close(&of->fh);
  Free(of->name);
  of->name = NULL;

  return iret == MPI_SUCCESS ? 0 : -1;
}

/***************************************************************/

static size_t
fwFileWriteOrdered(int fileID, const void *buffer, size_t len)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWO *of = openFiles + fileID;

  /* write buffer */
  xassert(len <= INT_MAX);
  MPI_Status status;
  xmpi(MPI_File_write_ordered(of->fh, (void *) buffer, (int) len, MPI_UNSIGNED_CHAR, &status));
  return len;
}

/***************************************************************/

static int
fcFileWriteOrdered(int fileID)
{
  assert(fileID >= 0 && (size_t) fileID < openFilesSize && openFiles[fileID].name);
  struct fileMPIFWO *of = openFiles + fileID;
  int iret = destroyAFiledataFileWriteOrdered(of);
  --openFilesFill;
  return iret;
}

/***************************************************************/

static int
fowFileWriteOrdered(const char *filename, const char *mode)
{
  if ((mode[0] != 'w' && mode[0] != 'W') || mode[0] == 0 || mode[1] != 0)
    xabort("Unsupported mode \"%s\" in parallel file open.", mode);

  struct cdiPioConf *conf = cdiPioGetConf();

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
  struct fileMPIFWO *of = openFiles + fileID;
  ++openFilesFill;
  initAFiledataFileWriteOrdered(of, filename, conf->writeAggBufLim);
  return (int) fileID;
}

/***************************************************************/

static void
finalizeFileWriteOrdered(void)
{
  if (openFilesFill)
    xabort("files still open on exit!");
  else
    Free(openFiles);
}

/***************************************************************/

void
cdiPioFileWriteOrderedInit(void)
{
  namespaceSwitchSet(NSSWITCH_FILE_OPEN, NSSW_FUNC(fowFileWriteOrdered));
  namespaceSwitchSet(NSSWITCH_FILE_CLOSE, NSSW_FUNC(fcFileWriteOrdered));
  namespaceSwitchSet(NSSWITCH_FILE_WRITE, NSSW_FUNC(fwFileWriteOrdered));
  namespaceSwitchSet(cdiPioExtraNSKeys[cdiPioEKFileWritingFinalize], NSSW_FUNC(finalizeFileWriteOrdered));
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
