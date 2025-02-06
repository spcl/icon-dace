#ifndef CDF_LAZY_GRID_H_
#define CDF_LAZY_GRID_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MMAP
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include <string.h>

#include "dmemory.h"
#include "cdf_int.h"
#include "grid.h"

struct xyValGet
{
  double scalefactor, addoffset;
  size_t start[3], count[3], size, dimsize;
  int datasetNCId, varNCId;
  short ndims;
};

struct cdfLazyGridIds
{
  int datasetNCId, varNCId;
};

struct cdfLazyGrid
{
  grid_t base;
  const struct gridVirtTable *baseVtable;
  struct cdfLazyGridIds cellAreaGet, xBoundsGet, yBoundsGet;
  struct xyValGet xValsGet, yValsGet;
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_t loadSerialize;
#endif
};

extern double *cdfPendingLoad;

void cdfLazyGridRenew(struct cdfLazyGrid *restrict *restrict gridpptr, int gridtype);
void cdfBaseGridRenew(struct cdfLazyGrid *restrict *restrict gridpptr, int gridtype);

void cdfLazyGridDestroy(struct cdfLazyGrid *lazyGrid);

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
