#ifndef CDI_PIO_DBUFFER_H
#define CDI_PIO_DBUFFER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

struct dBuffer
{
  size_t wr_pointer;
  size_t size;
  unsigned char *buffer;
};

void cdiPioDbufferInit(struct dBuffer *dbuffer, size_t bufSize);
bool cdiPioDbufferAppend(struct dBuffer *dbuffer, const void *data, size_t dataLen);
void cdiPioDbufferDestroy(struct dBuffer *dbuffer);

static inline int
cdiPioDbufferReset(struct dBuffer *dbuffer)
{
  dbuffer->wr_pointer = 0;
  return 0;
}

static inline size_t
cdiPioDbufferGetPos(struct dBuffer *dbuffer)
{
  return dbuffer->wr_pointer;
}

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
