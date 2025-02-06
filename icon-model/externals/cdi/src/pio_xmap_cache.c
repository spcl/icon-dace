#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include <mpi.h>
#include <yaxt.h>

#include "dmemory.h"
#include "error.h"
#include "pio_comm.h"
#include "pio_util.h"
#include "pio_xmap_cache.h"

enum
{
  LRAND48_BITS = 31,
};
#define MAX_CACHE_SIZE ((((size_t) 1) << LRAND48_BITS) - 1)

struct xmapCache
{
  size_t numUid, size, numEntries;
  Xt_uid *keys;
  int *numIndices;
  Xt_xmap xmaps[];
};

struct xmapCache *
cdiPioXmapCacheNew(size_t sizeEntries, size_t numUid)
{
  if (sizeEntries > MAX_CACHE_SIZE) xabort("cache cannot hold more than %zu entries", MAX_CACHE_SIZE);
  struct xmapCache *cache = Malloc(sizeof(struct xmapCache) + sizeof(Xt_xmap) * sizeEntries + sizeEntries);
  cache->numUid = numUid;
  cache->size = sizeEntries;
  cache->numEntries = 0;
  cache->keys = Malloc(sizeof(*cache->keys) * numUid * sizeEntries);
  cache->numIndices = Malloc(sizeof(*cache->numIndices) * numUid * sizeEntries);
  return cache;
}

size_t
cdiPioXmapCacheGetSize(struct xmapCache *cache)
{
  return cache ? cache->size : 0;
}

struct xmapCache *
cdiPioXmapCacheResize(struct xmapCache *cache, size_t sizeEntries)
{
  assert(cache);
  size_t cacheSize = cache->size;
  if (sizeEntries > MAX_CACHE_SIZE) xabort("cache cannot hold more than %zu entries", MAX_CACHE_SIZE);
  struct xmapCache *newCache;
  if (sizeEntries != cacheSize)
    {
      size_t numUid = cache->numUid;
      if (sizeEntries < cacheSize)
        {
          Xt_xmap *restrict xmaps = cache->xmaps;
          size_t numEntries = cache->numEntries;
          for (size_t i = sizeEntries; i < numEntries; ++i) xt_xmap_delete(xmaps[i]);
          cache->numEntries = numEntries > sizeEntries ? sizeEntries : numEntries;
        }
      cache->size = sizeEntries;
      newCache = Realloc(cache, sizeof(struct xmapCache) + sizeof(Xt_xmap) * sizeEntries + sizeEntries);
      size_t numKeysPerXmap = numUid * sizeEntries;
      newCache->keys = Realloc(newCache->keys, sizeof(*newCache->keys) * numKeysPerXmap);
      newCache->numIndices = Realloc(newCache->numIndices, sizeof(*newCache->numIndices) * numKeysPerXmap);
    }
  else
    newCache = cache;
  return newCache;
}

void
cdiPioXmapCacheDelete(struct xmapCache *cache)
{
  Xt_xmap *restrict xmaps = cache->xmaps;
  size_t numEntries = cache->numEntries;
  for (size_t i = 0; i < numEntries; ++i) xt_xmap_delete(xmaps[i]);
  Free(cache->numIndices);
  Free(cache->keys);
  Free(cache);
}

void
cdiPioXmapCacheAdd(struct xmapCache *cache, const Xt_uid *restrict keys, const int *restrict numIndices, Xt_xmap xmap)
{
  size_t cacheSize = cache->size, cacheFill = cache->numEntries, numUid = cache->numUid, cacheInsertPos;
  Xt_xmap *restrict xmaps = cache->xmaps;
  if (cacheFill < cacheSize)
    {
      cacheInsertPos = cacheFill;
      cache->numEntries = ++cacheFill;
    }
  else
    {
      /* we use lrand48 because it's guaranteed to generate 31 bits
       * per call */
      MPI_Comm collComm = commInqCommColl();
      long randVal = lrand48();
      xmpi(MPI_Bcast(&randVal, 1, MPI_LONG, 0, collComm));
      cacheInsertPos = (size_t) randVal % cacheSize;
      xt_xmap_delete(xmaps[cacheInsertPos]);
    }
  Xt_uid *restrict keys_ = cache->keys + numUid * cacheInsertPos;
  int *restrict numIndices_ = cache->numIndices + numUid * cacheInsertPos;
  for (size_t i = 0; i < numUid; ++i)
    {
      numIndices_[i] = numIndices[i];
      keys_[i] = keys[i];
    }
  xmaps[cacheInsertPos] = xmap;
}

Xt_xmap
cdiPioXmapCacheLookup(struct xmapCache *cache, const Xt_uid *restrict keys, int *restrict numIndices)
{
  size_t cacheSize = cache->size, cacheFill = cache->numEntries, numUid = cache->numUid;
  unsigned char *restrict matches = (unsigned char *) (cache->xmaps + cacheSize);
  for (size_t j = 0; j < cacheFill; ++j)
    {
      Xt_uid *restrict keys_ = cache->keys + j * numUid;
      Xt_uid match = 1;
      for (size_t i = 0; i < numUid; ++i) match &= (keys_[i] == keys[i]);
      matches[j] = (unsigned char) match;
      keys_ += numUid;
    }
  MPI_Comm collComm = commInqCommColl();
  xmpi(MPI_Allreduce(MPI_IN_PLACE, matches, (int) cacheFill, MPI_UNSIGNED_CHAR, MPI_LAND, collComm));
  for (size_t j = 0; j < cacheFill; ++j)
    if (matches[j])
      {
        Xt_xmap *restrict xmaps = cache->xmaps;
        int *restrict numIndices_ = cache->numIndices + j * numUid;
        for (size_t i = 0; i < numUid; ++i) numIndices[i] = numIndices_[i];
        return xmaps[j];
      }
  return NULL;
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
