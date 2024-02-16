#include "iterator_fallback.h"

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

struct CdiFallbackIterator
{
  CdiIterator super;
  char *path;  // needed for clone() & serialize()
  int streamId, vlistId, subtypeId;

  int variableCount, curVariable;
  int curLevelCount, curLevel;
  int curSubtypeCount, curSubtype;
  int curTimestep;
};

CdiIterator *
cdiFallbackIterator_getSuper(CdiFallbackIterator *me)
{
  return &me->super;
}

CdiIterator *
cdiFallbackIterator_new(const char *path, int filetype)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) Malloc(sizeof(*me));
  baseIterConstruct(&me->super, filetype);
  me->subtypeId = CDI_UNDEFID;  // Will be set in cdiFallbackIterator_nextField()
  me->curSubtypeCount = -1;     // Will be set in cdiFallbackIterator_nextField()
  me->curLevelCount = -1;       // Will be set in cdiFallbackIterator_nextField()
  // These values are chosen so that the natural increment at the start of cdiFallbackIterator_nextField() will correctly position
  // us at the first slice.
  me->curTimestep = 0;
  me->curVariable = -1;
  me->curSubtype = -1;
  me->curLevel = -1;
  me->streamId = streamOpenRead(path);
  if (me->streamId != CDI_UNDEFID)
    {
      me->vlistId = streamInqVlist(me->streamId);
      if (me->vlistId != CDI_UNDEFID && (me->variableCount = vlistNvars(me->vlistId)) > 0
          && streamInqTimestep(me->streamId, me->curTimestep) > 0 && (me->path = strdup(path)))
        {
          return (CdiIterator *) me;
        }
      Free(me->path);
      streamClose(me->streamId);
    }
  baseIterDestruct(&me->super);
  Free(me);
  return NULL;
}

void
cdiFallbackIterator_delete(CdiIterator *super)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  Free(me->path);
  streamClose(me->streamId);
  baseIterDestruct(super);
  Free(me);
}

// Fetches the info that is derived from the current variable. Most of this is published by the data members in the base class.
static void
fetchVariableInfo(CdiFallbackIterator *me)
{
  // Fetch data that's published via base class data members.
  me->super.datatype = vlistInqVarDatatype(me->vlistId, me->curVariable);
  me->super.timesteptype = vlistInqVarTsteptype(me->vlistId, me->curVariable);
  me->super.gridId = vlistInqVarGrid(me->vlistId, me->curVariable);
  int param = vlistInqVarParam(me->vlistId, me->curVariable);
  cdiDecodeParam(param, &me->super.param.number, &me->super.param.category, &me->super.param.discipline);

  // Fetch the current level and subtype counts.
  me->curLevelCount = zaxisInqSize(vlistInqVarZaxis(me->vlistId, me->curVariable));
  me->subtypeId = vlistInqVarSubtype(me->vlistId, me->curVariable);
  me->curSubtypeCount = (me->subtypeId == CDI_UNDEFID) ? 1 : subtypeInqSize(me->subtypeId);
}

CdiFallbackIterator *
cdiFallbackIterator_clone(CdiIterator *super)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;

  // Make another stream for this file. This yields an unadvanced iterator.
  CdiFallbackIterator *clone = (CdiFallbackIterator *) (void *) cdiFallbackIterator_new(me->path, me->super.filetype);
  if (clone)
    {
      // Point the clone to the same position in the file.
      clone->variableCount = me->variableCount;
      clone->curVariable = me->curVariable;
      clone->curLevelCount = me->curLevelCount;
      clone->curLevel = me->curLevel;
      clone->curSubtypeCount = me->curSubtypeCount;
      clone->curSubtype = me->curSubtype;
      clone->curTimestep = me->curTimestep;

      clone->super.isAdvanced = super->isAdvanced;
      if (super->isAdvanced) fetchVariableInfo(clone);
    }

  return clone;
}

char *
cdiFallbackIterator_serialize(CdiIterator *super)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;

  char *escapedPath = cdiEscapeSpaces(me->path);
  char *result = (char *) Malloc(strlen(escapedPath) + 7 * (3 * sizeof(int) * CHAR_BIT / 8 + 1) + 1);
  sprintf(result, "%s %d %d %d %d %d %d %d", escapedPath, me->variableCount, me->curVariable, me->curLevelCount, me->curLevel,
          me->curSubtypeCount, me->curSubtype, me->curTimestep);
  Free(escapedPath);
  return result;
}

CdiFallbackIterator *
cdiFallbackIterator_deserialize(const char *description)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) Malloc(sizeof(*me));
  if (!me) goto fail;

  description = baseIter_constructFromString(&me->super, description);

  while (*description == ' ') description++;
  me->path = cdiUnescapeSpaces(description, &description);
  if (!me->path) goto destructSuper;

  me->streamId = streamOpenRead(me->path);
  if (me->streamId == CDI_UNDEFID) goto freePath;
  me->vlistId = streamInqVlist(me->streamId);
  if (me->vlistId == CDI_UNDEFID) goto closeStream;

    // This reads one variable from the description string, does error checking, and advances the given string pointer.
#define decodeValue(variable, description)                                                                                      \
  do                                                                                                                            \
    {                                                                                                                           \
      const char *savedStart = description;                                                                                     \
      long long decodedValue                                                                                                    \
          = strtoll(description, (char **) &description, 0); /*The cast is a workaround for the wrong signature of strtoll().*/ \
      variable = (int) decodedValue;                                                                                            \
      if (savedStart == description) goto closeStream;                                                                          \
      if ((long long) decodedValue != (long long) variable) goto closeStream;                                                   \
    }                                                                                                                           \
  while (0)
  decodeValue(me->variableCount, description);
  decodeValue(me->curVariable, description);
  decodeValue(me->curLevelCount, description);
  decodeValue(me->curLevel, description);
  decodeValue(me->curSubtypeCount, description);
  decodeValue(me->curSubtype, description);
  decodeValue(me->curTimestep, description);
#undef decodeValue

  if (streamInqTimestep(me->streamId, me->curTimestep) <= 0) goto closeStream;
  if (me->super.isAdvanced) fetchVariableInfo(me);

  return me;

closeStream:
  streamClose(me->streamId);
freePath:
  Free(me->path);
destructSuper:
  baseIterDestruct(&me->super);
  Free(me);
fail:
  return NULL;
}

static int
advance(CdiFallbackIterator *me)
{
  me->curLevel++;
  if (me->curLevel >= me->curLevelCount)
    {
      me->curLevel = 0;
      me->curSubtype++;
      if (me->curSubtype >= me->curSubtypeCount)
        {
          me->curSubtype = 0;
          me->curVariable++;
          if (me->curVariable >= me->variableCount)
            {
              me->curVariable = 0;
              me->curTimestep++;
              if (streamInqTimestep(me->streamId, me->curTimestep) <= 0) return CDI_EEOF;
            }
        }
    }
  return CDI_NOERR;
}

int
cdiFallbackIterator_nextField(CdiIterator *super)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  int result = advance(me);
  if (result) return result;

  if (!me->curLevel && !me->curSubtype)
    fetchVariableInfo(me);  // Check whether we are processing a new variable/timestep and fetch the information that may have
                            // changed in this case.
  return CDI_NOERR;
}

char *
cdiFallbackIterator_inqTime(CdiIterator *super, CdiTimeType timeType)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;

  // retrieve the time information
  int taxisId = vlistInqTaxis(me->vlistId);
  int date = 0, time = 0;
  switch (timeType)
    {
    case kCdiTimeType_referenceTime:
      date = taxisInqRdate(taxisId);
      time = taxisInqRtime(taxisId);
      break;

    case kCdiTimeType_startTime:
      date = taxisInqVdate(taxisId);
      time = taxisInqVtime(taxisId);
      break;

    case kCdiTimeType_endTime:
      return NULL;  // The stream interface does not export the start/end times of statistical fields, so we treat all data as point
                    // of time data, returning the validity time as the start time.

    default: assert(0 && "internal error, please report this bug");
    }

  // decode the time information and reencode it into an ISO-compliant string
  int year, month, day, hour, minute, second;
  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);
  char *result = (char *) Malloc(4 + 1 + 2 + 1 + 2 + 1 + 2 + 1 + 2 + 1 + 2 + 4 + 1);
  sprintf(result, "%04d-%02d-%02dT%02d:%02d:%02d.000", year, month, day, hour, minute, second);
  return result;
}

int
cdiFallbackIterator_levelType(CdiIterator *super, int levelSelector, char **outName, char **outLongName, char **outStdName,
                              char **outUnit)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);
  (void) levelSelector;
#define copyString(outPointer, key)                                       \
  do                                                                      \
    {                                                                     \
      if (outPointer)                                                     \
        {                                                                 \
          char tempBuffer[CDI_MAX_NAME];                                  \
          int length = CDI_MAX_NAME;                                      \
          cdiInqKeyString(zaxisId, CDI_GLOBAL, key, tempBuffer, &length); \
          *outPointer = strdup(tempBuffer);                               \
        }                                                                 \
    }                                                                     \
  while (0)
  copyString(outName, CDI_KEY_NAME);
  copyString(outLongName, CDI_KEY_LONGNAME);
  copyString(outStdName, CDI_KEY_STDNAME);
  copyString(outUnit, CDI_KEY_UNITS);
#undef copyString
  int ltype = 0;
  cdiInqKeyInt(zaxisId, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
  return ltype;
}

int
cdiFallbackIterator_level(CdiIterator *super, int levelSelector, double *outValue1, double *outValue2)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);

  // handle NULL pointers once and for all
  double trash;
  if (!outValue1) outValue1 = &trash;
  if (!outValue2) outValue2 = &trash;

  // get the level value
  if (levelSelector)
    {
      *outValue1 = (zaxisInqLbounds(zaxisId, NULL)) ? zaxisInqLbound(zaxisId, me->curLevel) : zaxisInqLevel(zaxisId, me->curLevel);
    }
  else
    {
      *outValue1 = (zaxisInqUbounds(zaxisId, NULL)) ? zaxisInqUbound(zaxisId, me->curLevel) : zaxisInqLevel(zaxisId, me->curLevel);
    }
  *outValue2 = 0.0;

  // if this is a hybrid zaxis, lookup the coordinates in the vertical coordinate table
  ssize_t intLevel = (ssize_t) (2 * *outValue1);
  if (0 <= intLevel && intLevel < zaxisInqVctSize(zaxisId) - 1)
    {
      const double *coordinateTable = zaxisInqVctPtr(zaxisId);
      *outValue1 = coordinateTable[intLevel];
      *outValue2 = coordinateTable[intLevel + 1];
    }
  return CDI_NOERR;
}

int
cdiFallbackIterator_zaxisUuid(CdiIterator *super, int *outVgridNumber, int *outLevelCount, unsigned char outUuid[CDI_UUID_SIZE])
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  int zaxisId = vlistInqVarZaxis(me->vlistId, me->curVariable);
  int ltype = 0;
  cdiInqKeyInt(zaxisId, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
  if (ltype != ZAXIS_HYBRID) return CDI_EINVAL;
  if (outVgridNumber)
    {
      *outVgridNumber = 0;
      cdiInqKeyInt(zaxisId, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, outVgridNumber);
    }
  if (outLevelCount)
    {
      *outLevelCount = 0;
      cdiInqKeyInt(zaxisId, CDI_GLOBAL, CDI_KEY_NLEV, outLevelCount);
    }
  if (outUuid)
    {
      int length = CDI_UUID_SIZE;
      memset(outUuid, 0, length);
      cdiInqKeyBytes(zaxisId, CDI_GLOBAL, CDI_KEY_UUID, outUuid, &length);
    }
  return CDI_NOERR;
}

int
cdiFallbackIterator_inqTile(CdiIterator *super, int *outTileIndex, int *outTileAttribute)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
#ifndef __cplusplus
  if (!outTileIndex) outTileIndex = &(int){ 0 };
  if (!outTileAttribute) outTileAttribute = &(int){ 0 };
#else
  int dummy = 0;
  if (!outTileIndex) outTileIndex = &dummy;
  if (!outTileAttribute) outTileAttribute = &dummy;
#endif

  int error = CDI_NOERR;
  if (me->subtypeId == CDI_UNDEFID)  // must not call subtypeInqAttribute() with an invalid subtype ID, because it would abort the
                                     // program instead of returning an error
    {
      error = CDI_EINVAL;
    }
  else
    {
      if (subtypeInqAttribute(me->subtypeId, me->curSubtype, "tileIndex", outTileIndex)) error = CDI_EINVAL;
      if (subtypeInqAttribute(me->subtypeId, me->curSubtype, "tileAttribute", outTileAttribute)) error = CDI_EINVAL;
    }
  if (error) *outTileIndex = *outTileAttribute = -1;  // Guarantee defined values in case of an error.
  return error;
}

int
cdiFallbackIterator_inqTileCount(CdiIterator *super, int *outTileCount, int *outTileAttributeCount)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
#ifndef __cplusplus
  if (!outTileCount) outTileCount = &(int){ 0 };
  if (!outTileAttributeCount) outTileAttributeCount = &(int){ 0 };
#else
  int temp = 0;
  if (!outTileCount) outTileCount = &temp;
  if (!outTileAttributeCount) outTileAttributeCount = &temp;
#endif

  int error = CDI_NOERR;
  if (me->subtypeId == CDI_UNDEFID)  // must not call subtypeInqAttribute() with an invalid subtype ID, because it would abort the
                                     // program instead of returning an error
    {
      error = CDI_EINVAL;
    }
  else
    {
      if (subtypeInqAttribute(me->subtypeId, me->curSubtype, "numberOfTiles", outTileCount)) error = CDI_EINVAL;
      if (subtypeInqAttribute(me->subtypeId, me->curSubtype, "numberOfTileAttributes", outTileAttributeCount)) error = CDI_EINVAL;
    }
  if (error) *outTileCount = *outTileAttributeCount = -1;  // Guarantee defined values in case of an error.
  return CDI_NOERR;
}

char *
cdiFallbackIterator_copyVariableName(CdiIterator *super)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  return vlistCopyVarName(me->vlistId, me->curVariable);
}

void
cdiFallbackIterator_readField(CdiIterator *super, double *buffer, size_t *nmiss)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  SizeType missingValues = 0;
  streamReadVarSlice(me->streamId, me->curVariable, me->curLevel, buffer, &missingValues);
  if (nmiss) *nmiss = (size_t) missingValues;
}

void
cdiFallbackIterator_readFieldF(CdiIterator *super, float *buffer, size_t *nmiss)
{
  CdiFallbackIterator *me = (CdiFallbackIterator *) (void *) super;
  SizeType missingValues = 0;
  streamReadVarSliceF(me->streamId, me->curVariable, me->curLevel, buffer, &missingValues);
  if (nmiss) *nmiss = (size_t) missingValues;
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
