#ifndef CDI_PIO_XMAP_CACHE_H
#define CDI_PIO_XMAP_CACHE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include <yaxt.h>

struct xmapCache *cdiPioXmapCacheNew(size_t sizeEntries, size_t numUid);

size_t cdiPioXmapCacheGetSize(struct xmapCache *cache);

struct xmapCache *cdiPioXmapCacheResize(struct xmapCache *cache, size_t sizeEntries);

void cdiPioXmapCacheDelete(struct xmapCache *cache);

void cdiPioXmapCacheAdd(struct xmapCache *cache, const Xt_uid *restrict keys, const int *restrict numIndices, Xt_xmap xmap);

Xt_xmap cdiPioXmapCacheLookup(struct xmapCache *cache, const Xt_uid *restrict keys, int *restrict numIndices);

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
