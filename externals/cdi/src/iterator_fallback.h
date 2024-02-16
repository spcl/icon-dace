/*
 * A fallback implementation of the iterator interface that opens a stream under the hood.
 *
 * This implementation is mainly available to provide iterator access to file formats that don't support iterator access natively,
 * nevertheless, it allows the file to dictate the order in which data is read, possibly providing performance benefits.
 */

#ifndef INCLUDE_GUARD_CDI_ITERATOR_FALLBACK_H
#define INCLUDE_GUARD_CDI_ITERATOR_FALLBACK_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include "cdi.h"
#include "cdi_int.h"
#include "iterator.h"

typedef struct CdiFallbackIterator CdiFallbackIterator;

CdiIterator *cdiFallbackIterator_new(const char *path, int filetype);
CdiFallbackIterator *cdiFallbackIterator_clone(CdiIterator *me);
CdiIterator *cdiFallbackIterator_getSuper(CdiFallbackIterator *me);
char *cdiFallbackIterator_serialize(CdiIterator *me);
CdiFallbackIterator *cdiFallbackIterator_deserialize(const char *me);

int cdiFallbackIterator_nextField(CdiIterator *me);

char *cdiFallbackIterator_inqTime(CdiIterator *me, CdiTimeType timeType);
int cdiFallbackIterator_levelType(CdiIterator *me, int levelSelector, char **outName, char **outLongName, char **outStdName,
                                  char **outUnit);
int cdiFallbackIterator_level(CdiIterator *me, int levelSelector, double *outValue1, double *outValue2);
int cdiFallbackIterator_zaxisUuid(CdiIterator *me, int *outVgridNumber, int *outLevelCount, unsigned char outUuid[CDI_UUID_SIZE]);
char *cdiFallbackIterator_copyVariableName(CdiIterator *me);
int cdiFallbackIterator_inqTile(CdiIterator *me, int *outTileIndex, int *outTileAttribute);
int cdiFallbackIterator_inqTileCount(CdiIterator *me, int *outTileCount, int *outTileAttributeCount);

void cdiFallbackIterator_readField(CdiIterator *me, double *buffer, size_t *nmiss);
void cdiFallbackIterator_readFieldF(CdiIterator *me, float *buffer, size_t *nmiss);

void cdiFallbackIterator_delete(CdiIterator *super);

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
