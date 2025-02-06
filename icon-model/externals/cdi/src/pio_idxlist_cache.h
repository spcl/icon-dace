#ifndef CDI_PIO_IDXLIST_CACHE_H
#define CDI_PIO_IDXLIST_CACHE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

struct cdiPioIdxlistCache *cdiPioIdxlistCacheNew(size_t sizeEntries);

size_t cdiPioIdxlistCacheGetSize(struct cdiPioIdxlistCache *cache);

struct cdiPioIdxlistCache *cdiPioIdxlistCacheResize(struct cdiPioIdxlistCache *cache, size_t sizeEntries);

void cdiPioIdxlistCacheDelete(struct cdiPioIdxlistCache *cache);

Xt_idxlist cdiPioIdxlistCacheAddSection2D(struct cdiPioIdxlistCache *cache, const Xt_int wholeShape[2], const Xt_int sliceOrigin[2],
                                          const int sliceShape[2]);

Xt_idxlist cdiPioIdxlistCacheAddSection3D(struct cdiPioIdxlistCache *cache, const Xt_int wholeShape[3], const Xt_int sliceOrigin[3],
                                          const int sliceShape[3]);

Xt_idxlist cdiPioIdxlistCacheAddStripes1(struct cdiPioIdxlistCache *cache, Xt_int start, int nstrides);

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
