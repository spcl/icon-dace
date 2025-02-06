#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "dmemory.h"

#include "pio_dbuffer.h"
#include "pio_util.h"

void
cdiPioDbufferInit(struct dBuffer *db, size_t bufSize)
{
  assert(db);
  db->wr_pointer = 0;

#ifndef _SX
  size_t pagesize = (size_t) (sysconf(_SC_PAGESIZE));

  xdebug("cdiPioDbufferInit(): pagesize = %zu bytes, size = %zu", pagesize, bufSize);

  db->size = bufSize = (bufSize + pagesize - 1) / pagesize * pagesize;
  void *buf = NULL;
  int status;
  if ((status = posix_memalign(&buf, pagesize, bufSize)) != 0)
    {
      switch (status)
        {
        case EINVAL:
          xabort("The alignment argument was not a power of two,"
                 " or was not a multiple of sizeof(void *).");
#ifndef __GNUC__
          break;
#endif
        case ENOMEM:
          xabort("There was insufficient memory to fulfill the"
                 " allocation request.");
#ifndef __GNUC__
          break;
#endif
        }
    }
  db->buffer = (unsigned char *) buf;
#else
  db->size = bufSize;
  db->buffer = (unsigned char *) Malloc(bufSize);
#endif
}

void
cdiPioDbufferDestroy(struct dBuffer *db)
{
  free(db->buffer);
}

static size_t
dbufferFreesize(struct dBuffer *dbuffer)
{
  size_t free_size = (size_t) (dbuffer->size - cdiPioDbufferGetPos(dbuffer));
  return free_size;
}

bool
cdiPioDbufferAppend(struct dBuffer *dbuffer, const void *data, size_t dataLen)
{
  size_t space_left = dbufferFreesize(dbuffer);
  bool enoughSpaceLeft = dataLen <= space_left;
  if (dataLen <= space_left)
    {
      size_t wr_ptr = dbuffer->wr_pointer;
      memcpy(dbuffer->buffer + wr_ptr, data, dataLen);
      dbuffer->wr_pointer = wr_ptr + dataLen;
    }
  return !enoughSpaceLeft;
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
