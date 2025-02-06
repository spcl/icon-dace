#ifndef PIO_ID_LIST_H
#define PIO_ID_LIST_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdlib.h>

#include "cdi.h"
#include "error.h"
#include "dmemory.h"

struct idList
{
  int *entries;
  size_t size;
};

#define emptyIDList ((struct idList){ NULL, 0 })

static inline size_t
indexOfID(const struct idList *list, int ID)
{
  size_t index = SIZE_MAX;
  for (size_t i = 0; i < list->size; ++i)
    if (list->entries[i] == ID) index = i;
  return index;
}

static inline size_t
insertID(struct idList *list, int ID)
{
  size_t index = indexOfID(list, CDI_UNDEFID);
  if (index == SIZE_MAX)
    {
      index = list->size;
      list->entries = Realloc(list->entries, ++list->size * sizeof(list->entries[0]));
    }
  list->entries[index] = ID;
  return index;
}

static inline void
removeID(struct idList *list, int ID)
{
  size_t index = indexOfID(list, ID);
  xassert(index != SIZE_MAX);
  list->entries[index] = CDI_UNDEFID;
}

static inline void
idSetDestroy(struct idList *list)
{
  free(list->entries);
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
