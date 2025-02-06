#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdlib.h>

#include <yaxt.h>

#include "dmemory.h"
#include "error.h"
#include "pio_idxlist_cache.h"

enum
{
  LRAND48_BITS = 31,
};
#define MAX_CACHE_SIZE ((((size_t) 1) << LRAND48_BITS) - 1)

enum cdiPioIdxlistType
{
  CDIPIO_IDXLISTTYPE_EMPTY,
  CDIPIO_IDXLISTTYPE_SECTION2D,
  CDIPIO_IDXLISTTYPE_SECTION3D,
  CDIPIO_IDXLISTTYPE_STRIPES1,
};

struct sectionNDDesc
{
  Xt_int wholeShape[3], sliceOrigin[3];
  int sliceShape[3];
};

struct stripes1Desc
{
  Xt_int start;
  int nstrides;
};

struct idxlistDesc
{
  enum cdiPioIdxlistType type;
  union
  {
    struct sectionNDDesc sectionND;
    struct stripes1Desc stripes1;
  } desc;
  Xt_idxlist idxlist;
};

struct cdiPioIdxlistCache
{
  size_t size, numEntries;
  struct idxlistDesc entries[];
};

struct cdiPioIdxlistCache *
cdiPioIdxlistCacheNew(size_t sizeEntries)
{
  if (sizeEntries > MAX_CACHE_SIZE) xabort("cache cannot hold more than %zu entries", MAX_CACHE_SIZE);
  struct cdiPioIdxlistCache *cache = Malloc(sizeof(struct cdiPioIdxlistCache) + sizeof(struct idxlistDesc) * sizeEntries);
  cache->size = sizeEntries;
  cache->numEntries = 0;
  struct idxlistDesc *restrict entries = cache->entries;
  for (size_t i = 0; i < sizeEntries; ++i) entries[i].type = CDIPIO_IDXLISTTYPE_EMPTY;
  return cache;
}

size_t
cdiPioIdxlistCacheGetSize(struct cdiPioIdxlistCache *cache)
{
  return cache ? cache->size : 0;
}

struct cdiPioIdxlistCache *
cdiPioIdxlistCacheResize(struct cdiPioIdxlistCache *cache, size_t sizeEntries)
{
  size_t cacheSize = cache ? cache->size : 0;
  if (sizeEntries > MAX_CACHE_SIZE) xabort("cache cannot hold more than %zu entries", MAX_CACHE_SIZE);
  struct cdiPioIdxlistCache *newCache;
  if (sizeEntries < cacheSize)
    {
      struct idxlistDesc *restrict entries = cache->entries;
      size_t numEntries = cache->numEntries;
      for (size_t i = sizeEntries; i < numEntries; ++i) xt_idxlist_delete(entries[i].idxlist);
      cache->numEntries = numEntries > sizeEntries ? sizeEntries : numEntries;
      cache->size = sizeEntries;
      newCache = Realloc(cache, (sizeEntries ? sizeof(struct cdiPioIdxlistCache) : 0) + sizeof(struct idxlistDesc) * sizeEntries);
    }
  else if (sizeEntries > cacheSize)
    {
      newCache = Realloc(cache, sizeof(struct cdiPioIdxlistCache) + sizeof(struct idxlistDesc) * sizeEntries);
      if (!cache) newCache->numEntries = 0;
      newCache->size = sizeEntries;
      struct idxlistDesc *restrict entries = newCache->entries;
      for (size_t i = cacheSize; i < sizeEntries; ++i) entries[i].type = CDIPIO_IDXLISTTYPE_EMPTY;
    }
  else
    newCache = cache;
  return newCache;
}

void
cdiPioIdxlistCacheDelete(struct cdiPioIdxlistCache *cache)
{
  if (cache)
    {
      size_t numEntries = cache->numEntries;
      struct idxlistDesc *restrict entries = cache->entries;
      for (size_t i = 0; i < numEntries; ++i) xt_idxlist_delete(entries[i].idxlist);
      Free(cache);
    }
}

static Xt_idxlist
cdiPioIdxlistCacheAddSectionND(struct cdiPioIdxlistCache *cache, enum cdiPioIdxlistType type, size_t ndims,
                               const Xt_int wholeShape[], const Xt_int sliceOrigin[], const int sliceShape[])
{
  size_t cacheSize = cache->size, cacheFill = cache->numEntries, cacheInsertPos;
  struct idxlistDesc *restrict entries = cache->entries;
  for (size_t i = 0; i < cacheFill; ++i)
    if (entries[i].type == type)
      {
        bool equal = true;
        for (size_t j = 0; j < ndims; ++j)
          equal &= ((entries[i].desc.sectionND.wholeShape[j] == wholeShape[j])
                    & (entries[i].desc.sectionND.sliceOrigin[j] == sliceOrigin[j])
                    & (entries[i].desc.sectionND.sliceShape[j] == sliceShape[j]));
        if (equal) return entries[i].idxlist;
      }
  /* at this point it's definite that the requested section is not
   * in the cache, determine to replace or append is next */
  if (cacheFill < cacheSize)
    {
      cacheInsertPos = cacheFill;
      cache->numEntries = ++cacheFill;
    }
  else
    {
      /* we use lrand48 because it's guaranteed to generate 31 bits
       * per call */
      cacheInsertPos = (size_t) (lrand48());
      cacheInsertPos %= cacheSize;
      xt_idxlist_delete(entries[cacheInsertPos].idxlist);
    }
  entries[cacheInsertPos].type = type;
  for (size_t j = 0; j < ndims; ++j)
    {
      entries[cacheInsertPos].desc.sectionND.wholeShape[j] = wholeShape[j];
      entries[cacheInsertPos].desc.sectionND.sliceOrigin[j] = sliceOrigin[j];
      entries[cacheInsertPos].desc.sectionND.sliceShape[j] = sliceShape[j];
    }
  entries[cacheInsertPos].idxlist = xt_idxsection_new(0, (int) ndims, wholeShape, sliceShape, sliceOrigin);
  return entries[cacheInsertPos].idxlist;
}

Xt_idxlist
cdiPioIdxlistCacheAddSection2D(struct cdiPioIdxlistCache *cache, const Xt_int wholeShape[2], const Xt_int sliceOrigin[2],
                               const int sliceShape[2])
{
  return cdiPioIdxlistCacheAddSectionND(cache, CDIPIO_IDXLISTTYPE_SECTION2D, 2, wholeShape, sliceOrigin, sliceShape);
}

Xt_idxlist
cdiPioIdxlistCacheAddSection3D(struct cdiPioIdxlistCache *cache, const Xt_int wholeShape[3], const Xt_int sliceOrigin[3],
                               const int sliceShape[3])
{
  return cdiPioIdxlistCacheAddSectionND(cache, CDIPIO_IDXLISTTYPE_SECTION3D, 3, wholeShape, sliceOrigin, sliceShape);
}

Xt_idxlist
cdiPioIdxlistCacheAddStripes1(struct cdiPioIdxlistCache *cache, Xt_int start, int nstrides)
{
  size_t cacheSize = cache->size, cacheFill = cache->numEntries, cacheInsertPos;
  struct idxlistDesc *restrict entries = cache->entries;
  for (size_t i = 0; i < cacheFill; ++i)
    if (entries[i].type == CDIPIO_IDXLISTTYPE_STRIPES1 && entries[i].desc.stripes1.start == start
        && entries[i].desc.stripes1.nstrides == nstrides)
      return entries[i].idxlist;
  /* at this point it's definite that the requested section is not
   * in the cache, determine to replace or append is next */
  if (cacheFill < cacheSize)
    {
      cacheInsertPos = cacheFill;
      cache->numEntries = ++cacheFill;
    }
  else
    {
      /* we use lrand48 because it's guaranteed to generate 31 bits
       * per call */
      cacheInsertPos = (size_t) (lrand48());
      cacheInsertPos %= cacheSize;
      xt_idxlist_delete(entries[cacheInsertPos].idxlist);
    }
  entries[cacheInsertPos].type = CDIPIO_IDXLISTTYPE_STRIPES1;
  entries[cacheInsertPos].desc.stripes1.start = start;
  entries[cacheInsertPos].desc.stripes1.nstrides = nstrides;
  struct Xt_stripe stripe = { .start = start, .stride = 1, .nstrides = nstrides };
  entries[cacheInsertPos].idxlist = xt_idxstripes_new(&stripe, 1);
  return entries[cacheInsertPos].idxlist;
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
