/*
 * An implementation of the iterator interface for GRIB files.
 * Since GRIB files do not contain an index, this avoids scanning the entire file to generate an in-memory index as streamOpenRead()
 * does. Consequently, using this interface is much more efficient for GRIB files.
 */

#ifndef INCLUDE_GUARD_CDI_ITERATOR_GRIB_H
#define INCLUDE_GUARD_CDI_ITERATOR_GRIB_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "iterator.h"
#include "input_file.h"

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h>
#endif

typedef struct recordList recordList;

CdiIterator *cdiGribIterator_new(const char *path, int filetype);
CdiGribIterator *cdiGribIterator_makeClone(CdiIterator *me);
CdiIterator *cdiGribIterator_getSuper(CdiGribIterator *me);
char *cdiGribIterator_serialize(CdiIterator *me);
CdiGribIterator *cdiGribIterator_deserialize(const char *me);

int cdiGribIterator_nextField(CdiIterator *me);

char *cdiGribIterator_inqTime(CdiIterator *me, CdiTimeType timeType);
int cdiGribIterator_levelType(CdiIterator *me, int levelSelector, char **outName, char **outLongName, char **outStdName,
                              char **outUnit);
int cdiGribIterator_level(CdiIterator *me, int levelSelector, double *outValue1, double *outValue2);
int cdiGribIterator_zaxisUuid(CdiIterator *me, int *outVgridNumber, int *outLevelCount, unsigned char outUuid[CDI_UUID_SIZE]);
int cdiGribIterator_inqTile(CdiIterator *me, int *outTileIndex, int *outTileAttribute);
int cdiGribIterator_inqTileCount(CdiIterator *me, int *outTileCount, int *outTileAttributeCount);
char *cdiGribIterator_copyVariableName(CdiIterator *me);

void cdiGribIterator_readField(CdiIterator *me, double *buffer, size_t *numMissVals);
void cdiGribIterator_readFieldF(CdiIterator *me, float *buffer, size_t *numMissVals);

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
