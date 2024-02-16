#ifndef INCLUDE_GUARD_CDI_GRIBAPI_UTILITIES_H
#define INCLUDE_GUARD_CDI_GRIBAPI_UTILITIES_H

#ifdef HAVE_LIBGRIB_API

#include "cdi_int.h"
#include "grid.h"

#include <grib_api.h>

#include <stdbool.h>

char *gribCopyString(grib_handle *gribHandle, const char *key);
bool gribCheckString(grib_handle *gribHandle, const char *key, const char *expectedValue);

bool gribCheckLong(grib_handle *gribHandle, const char *key, long expectedValue);
long gribGetLong(grib_handle *gh, const char *key);
long gribGetLongDefault(grib_handle *gribHandle, const char *key, long defaultValue);

double gribGetDouble(grib_handle *gh, const char *key);
double gribGetDoubleDefault(grib_handle *gribHandle, const char *key, double defaultValue);

size_t gribGetArraySize(grib_handle *gribHandle, const char *key);
void gribGetDoubleArray(grib_handle *gribHandle, const char *key,
                        double *array);  // The caller is responsible to ensure a sufficiently large buffer.
void gribGetLongArray(grib_handle *gribHandle, const char *key,
                      long *array);  // The caller is responsible to ensure a sufficiently large buffer.

long gribEditionNumber(grib_handle *gh);
char *gribMakeTimeString(grib_handle *gh, CdiTimeType timeType);  // Returns NULL if timeType is kCdiTimeType_endTime and the field
                                                                  // does not have an integration period (statistical data).
int gribapiTimeIsFC(grib_handle *gh);
int gribapiGetTsteptype(grib_handle *gh);
int gribGetDatatype(grib_handle *gribHandle);
int gribapiGetParam(grib_handle *gh);
int gribapiGetGridType(grib_handle *gh);
bool gribapiGetGrid(grib_handle *gh, grid_t *grid);
size_t gribapiGetGridsize(grib_handle *gh);

#ifdef HIRLAM_EXTENSIONS
void gribapiSetDataTimeRangeIndicator(grib_handle *gh, int timeRangeIndicator);
void gribapiGetDataTimeRangeIndicator(grib_handle *gh, int *timeRangeIndicator);
#endif  // #ifdef HIRLAM_EXTENSIONS

extern struct cdiGribAPI_ts_str_map_elem
{
  long productionTemplate;
  const char sname[8];
} cdiGribAPI_ts_str_map[];

#endif

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
