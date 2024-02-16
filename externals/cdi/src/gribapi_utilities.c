#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBGRIB_API

#include "gribapi_utilities.h"

#include "cdi.h"
#include "dmemory.h"
#include "error.h"
#include "gribapi.h"
#include "grid.h"

#include <assert.h>
#include <time.h>

#define FAIL_ON_GRIB_ERROR(function, gribHandle, key, ...)                                                                 \
  do                                                                                                                       \
    {                                                                                                                      \
      int errorCode = (int) function(gribHandle, key, __VA_ARGS__);                                                        \
      if (errorCode)                                                                                                       \
        {                                                                                                                  \
          fprintf(stderr, "%s:%d: Error in function `%s`: `%s` returned error code %d for key \"%s\"", __FILE__, __LINE__, \
                  __func__, #function, errorCode, key);                                                                    \
          exit(errorCode);                                                                                                 \
        }                                                                                                                  \
    }                                                                                                                      \
  while (0)

// A simple wrapper for grib_get_string() that returns a newly allocated string.
char *
gribCopyString(grib_handle *gribHandle, const char *key)
{
  size_t length;
#ifdef HAVE_GRIB_GET_LENGTH
  if (!grib_get_length(gribHandle, key, &length))
    {
      char *result = (char *) Malloc(length);
      if (!grib_get_string(gribHandle, key, result, &length))
        result = (char *) Realloc(result, length);
      else
        {
          Free(result);
          result = NULL;
        }
      return result;
    }
  else
    return NULL;
#else
  length = 1024; /* there's an implementation limit
                  * that makes strings longer than
                  * this unlikely in grib_api versions
                  * not providing grib_get_length */
  int rc;
  char *result = (char *) Malloc(length);
  while ((rc = grib_get_string(gribHandle, key, result, &length)) == GRIB_BUFFER_TOO_SMALL || rc == GRIB_ARRAY_TOO_SMALL)
    {
      if (length <= 1024UL * 1024UL)
        {
          length *= 2;
          result = Realloc(result, length);
        }
      else
        break;
    }
  if (!rc)
    result = Realloc(result, length);
  else
    {
      Free(result);
      result = NULL;
    }
  return result;
#endif
}

// A simple wrapper for grib_get_string() for the usecase that the result is only compared to a given constant string.
// Returns true if the key exists and the value is equal to the given string.
bool
gribCheckString(grib_handle *gribHandle, const char *key, const char *expectedValue)
{
  size_t expectedLength = strlen(expectedValue) + 1;
#ifdef HAVE_GRIB_GET_LENGTH
  size_t length;
  if (grib_get_length(gribHandle, key, &length)) return false;
  if (length != expectedLength) return false;
  char *value = (char *) Malloc(length);
  if (grib_get_string(gribHandle, key, value, &length)) return false;
  int rc = str_is_equal(value, expectedValue);
  Free(value);
#else
  char *value = gribCopyString(gribHandle, key);
  int rc = value ? (strlen(value) + 1 == expectedLength ? str_is_equal(value, expectedValue) : false) : false;
  Free(value);
#endif
  return rc;
}

// A simple wrapper for grib_get_long() for the usecase that the result is only compared to a given constant value.
// Returns true if the key exists and the value is equal to the given one.
bool
gribCheckLong(grib_handle *gribHandle, const char *key, long expectedValue)
{
  long value;
  if (grib_get_long(gribHandle, key, &value)) return false;
  return value == expectedValue;
}

// A simple wrapper for grib_get_long() for the usecase that failure to fetch the value is fatal.
long
gribGetLong(grib_handle *gh, const char *key)
{
  long result;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, key, &result);
  return result;
}

// A simple wrapper for grib_get_long() for the usecase that a default value is used in the case that the operation fails.
long
gribGetLongDefault(grib_handle *gribHandle, const char *key, long defaultValue)
{
  long result;
  if (grib_get_long(gribHandle, key, &result) || result == GRIB_MISSING_LONG) result = defaultValue;
  return result;
}

// A simple wrapper for grib_get_double() for the usecase that failure to fetch the value is fatal.
double
gribGetDouble(grib_handle *gh, const char *key)
{
  double result;
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, key, &result);
  return result;
}

// A sample wrapper for grib_get_double() for the usecase that a default value is used in the case that the operation fails.
double
gribGetDoubleDefault(grib_handle *gribHandle, const char *key, double defaultValue)
{
  double result;
  if (grib_get_double(gribHandle, key, &result) || IS_EQUAL(result, GRIB_MISSING_DOUBLE)) result = defaultValue;
  return result;
}

// A simple wrapper for grib_get_size() for the usecase that failure to fetch the value is fatal.
size_t
gribGetArraySize(grib_handle *gribHandle, const char *key)
{
  size_t result;
  FAIL_ON_GRIB_ERROR(grib_get_size, gribHandle, key, &result);
  return result;
}

// A simple wrapper for grib_get_double_array() for the usecase that failure to fetch the data is fatal.
void
gribGetDoubleArray(grib_handle *gribHandle, const char *key, double *array)
{
  size_t valueCount = gribGetArraySize(gribHandle, key);
  FAIL_ON_GRIB_ERROR(grib_get_double_array, gribHandle, key, array, &valueCount);
}

// A simple wrapper for grib_get_long_array() for the usecase that failure to fetch the data is fatal.
void
gribGetLongArray(grib_handle *gribHandle, const char *key, long *array)
{
  size_t valueCount = gribGetArraySize(gribHandle, key);
  FAIL_ON_GRIB_ERROR(grib_get_long_array, gribHandle, key, array, &valueCount);
}

// We need the edition number so frequently, that it's convenient to give it its own function.
long
gribEditionNumber(grib_handle *gh)
{
  return gribGetLong(gh, "editionNumber");
}

// This return value of this should be passed to a call to resetTz(), it is a malloc'ed string with the content of the TZ
// environment variable before the call (or NULL if that was not set).
static char *
setUtc(void)
{
  char *temp = getenv("TZ"), *result = NULL;
  if (temp) result = strdup(temp);
  setenv("TZ", "UTC", 1);
  return result;
}

// Undoes the effect of setUtc(), pass to it the return value of the corresponding setUtc() call, it will free the string.
static void
resetTz(char *savedTz)
{
  if (savedTz)
    {
      setenv("TZ", savedTz, 1);
      Free(savedTz);
    }
  else
    {
      unsetenv("TZ");
    }
}

// This function uses the system functions to normalize the date representation according to the gregorian calendar.
// Returns zero on success.
static int
normalizeDays(struct tm *me)
{
  char *savedTz = setUtc();  // Ensure that mktime() does not interprete the date according to our local time zone.

  int result = (mktime(me) == (time_t) -1);  // This does all the heavy lifting.

  resetTz(savedTz);
  return result;
}

// Returns zero on success.
static int
addSecondsToDate(struct tm *me, long long amount)
{
  // It is irrelevant here whether days are zero or one based, the correction would have be undone again so that it is effectless.
  long long seconds = ((me->tm_mday * 24ll + me->tm_hour) * 60 + me->tm_min) * 60
                      + me->tm_sec;  // The portion of the date that uses fixed increments.
  seconds += amount;
  me->tm_mday = (int) (seconds / 24 / 60 / 60);
  seconds -= (long long) me->tm_mday * 24 * 60 * 60;
  me->tm_hour = (int) (seconds / 60 / 60);
  seconds -= (long long) me->tm_hour * 60 * 60;
  me->tm_min = (int) (seconds / 60);
  seconds -= (long long) (me->tm_min * 60);
  me->tm_sec = (int) seconds;
  return normalizeDays(me);
}

static void
addMonthsToDate(struct tm *me, long long amount)
{
  long long months = me->tm_year * 12ll + me->tm_mon;
  months += amount;
  me->tm_year = (int) (months / 12);
  months -= (long long) me->tm_year * 12;
  me->tm_mon = (int) months;
}

// unit is a value according to code table 4.4 of the GRIB2 specification, returns non-zero on error
static int
addToDate(struct tm *me, long long amount, long unit)
{
  switch (unit)
    {
    case 0: return addSecondsToDate(me, 60 * amount);            // minute
    case 1: return addSecondsToDate(me, 60 * 60 * amount);       // hour
    case 2: return addSecondsToDate(me, 24 * 60 * 60 * amount);  // day

    case 3: addMonthsToDate(me, amount); return 0;             // month
    case 4: addMonthsToDate(me, 12 * amount); return 0;        // year
    case 5: addMonthsToDate(me, 10 * 12 * amount); return 0;   // decade
    case 6: addMonthsToDate(me, 30 * 12 * amount); return 0;   // normal
    case 7: addMonthsToDate(me, 100 * 12 * amount); return 0;  // century

    case 10: return addSecondsToDate(me, 3 * 60 * 60 * amount);   // eighth of a day
    case 11: return addSecondsToDate(me, 6 * 60 * 60 * amount);   // quarter day
    case 12: return addSecondsToDate(me, 12 * 60 * 60 * amount);  // half day
    case 13: return addSecondsToDate(me, amount);                 // second

    default: return 1;  // reserved, unknown, or missing
    }
}

static char *
makeDateString(struct tm *me)
{
  const size_t length = 4 + 1 + 2 + 1 + 2 + 1 + 2 + 1 + 2 + 1 + 2 + 4 + 1;
  char *result = (char *) Malloc(length);
  snprintf(result, length, "%04d-%02d-%02dT%02d:%02d:%02d.000", me->tm_year + 1900, me->tm_mon + 1, me->tm_mday, me->tm_hour,
           me->tm_min, me->tm_sec);
  return result;
}

// FIXME: This ignores any calendar definition that might be present.
// XXX: Identification templates are not implemented in grib_api-1.12.3, so even if I implemented the other calendars now, it
// wouldn't be possible to use them.
static int
getAvailabilityOfRelativeTimes(grib_handle *gh, bool *outHaveForecastTime, bool *outHaveTimeRange)
{
  switch (gribGetLong(gh, "productDefinitionTemplateNumber"))
    {
    case 20:
    case 30:
    case 31:
    case 254:
    case 311:
    case 2000: *outHaveForecastTime = false, *outHaveTimeRange = false; return 0;

    // case 55 and case 40455 are the same: 55 is the proposed standard value, 40455 is the value in the local use range that is
    // used by the dwd until the standard is updated.
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 15:
    case 32:
    case 33:
    case 40:
    case 41:
    case 44:
    case 45:
    case 48:
    case 51:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:
    case 58:
    case 60:
    case 1000:
    case 1002:
    case 1100:
    case 40033:
    case 40455:
    case 40456: *outHaveForecastTime = true, *outHaveTimeRange = false; return 0;

    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 34:
    case 42:
    case 43:
    case 46:
    case 47:
    case 61:
    case 67:
    case 68:
    case 91:
    case 1001:
    case 1101:
    case 40034: *outHaveForecastTime = true, *outHaveTimeRange = true; return 0;

    default: return 1;
    }
}

char *
gribMakeTimeString(grib_handle *gh, CdiTimeType timeType)
{
  // Get the parts of the reference date.
  struct tm date;
  date.tm_mon = (int) gribGetLong(gh, "month") - 1;  // months are zero based in struct tm and one based in GRIB
  date.tm_mday = (int) gribGetLong(gh, "day");
  date.tm_hour = (int) gribGetLong(gh, "hour");
  date.tm_min = (int) gribGetLong(gh, "minute");
  date.tm_isdst = 0;

  if (gribEditionNumber(gh) == 1)
    {
      date.tm_year = (int) gribGetLong(gh, "yearOfCentury");  // years are -1900 based both in struct tm and GRIB1
    }
  else
    {
      date.tm_year = (int) gribGetLong(gh, "year") - 1900;  // years are -1900 based in struct tm and zero based in GRIB2
      date.tm_sec = (int) gribGetLong(gh, "second");

      // If the start or end time are requested, we need to take the relative times into account.
      if (timeType != kCdiTimeType_referenceTime)
        {
          // Determine whether we have a forecast time and a time range.
          bool haveForecastTime, haveTimeRange;
          if (getAvailabilityOfRelativeTimes(gh, &haveForecastTime, &haveTimeRange)) return NULL;
          if (timeType == kCdiTimeType_endTime && !haveTimeRange)
            return NULL;  // tell the caller that the requested time does not exist

          // If we have relative times, apply the relative times to the date
          if (haveForecastTime)
            {
              long offset = gribGetLongDefault(gh, "forecastTime", 0);
              // if (stepUnits == indicatorOfUnitOfTimeRange) assert(startStep == forecastTime)
              long offsetUnit = gribGetLongDefault(gh, "indicatorOfUnitOfTimeRange", 255);
              if (addToDate(&date, offset, offsetUnit)) return NULL;
              if (timeType == kCdiTimeType_endTime)
                {
                  assert(haveTimeRange);
                  long range = gribGetLongDefault(gh, "lengthOfTimeRange", 0);
                  // if (stepUnits == indicatorOfUnitForTimeRange) assert(endStep == startStep + lengthOfTimeRange)
                  long rangeUnit = gribGetLongDefault(gh, "indicatorOfUnitForTimeRange", 255);
                  if (addToDate(&date, range, rangeUnit)) return NULL;
                }
            }
        }
    }

  // Bake the date into a string.
  return makeDateString(&date);
}

int
gribapiTimeIsFC(grib_handle *gh)
{
  if (gribEditionNumber(gh) <= 1) return true;

  long sigofrtime;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "significanceOfReferenceTime", &sigofrtime);
  return sigofrtime != 3;
}

struct cdiGribAPI_ts_str_map_elem cdiGribAPI_ts_str_map[] = {
  // clang-format off
  [TSTEP_INSTANT] = {  0, "instant" },
  [TSTEP_AVG]     = {  8, "avg" },
  [TSTEP_ACCUM]   = {  8, "accum" },
  [TSTEP_MAX]     = {  8, "max" },
  [TSTEP_MIN]     = {  8, "min" },
  [TSTEP_DIFF]    = {  8, "diff" },
  [TSTEP_RMS]     = {  8, "rms" },
  [TSTEP_SD]      = {  8, "sd" },
  [TSTEP_COV]     = {  8, "cov" },
  [TSTEP_RATIO]   = {  8, "ratio" },
  [TSTEP_SUM]     = {  8, "sum" },
                    {  0, "" }
  // clang-format on
};

// Fetches the value of the "stepType" key and converts it into a constant in the TSTEP_* range.
int
gribapiGetTsteptype(grib_handle *gh)
{
  size_t len = 256;
  char stepType[256];
  int tsteptype = TSTEP_INSTANT;
  static bool lprint = true;

  if (gribapiTimeIsFC(gh))
    {
      int status = grib_get_string(gh, "stepType", stepType, &len);
      if (status == 0 && len > 1 && len < 256)
        {
          for (int i = TSTEP_INSTANT; cdiGribAPI_ts_str_map[i].sname[0]; ++i)
            if (strncmp(cdiGribAPI_ts_str_map[i].sname, stepType, len) == 0)
              {
                tsteptype = i;
                goto tsteptypeFound;
              }

          if (lprint)
            {
              Message("Time stepType %s unsupported, set to instant!", stepType);
              lprint = false;
            }
          // printf("stepType: %s %ld %d\n", stepType, len, tsteptype);
        }

      long typeOfStat;
      status = grib_get_long(gh, "typeOfStatisticalProcessing", &typeOfStat);
      if (status == 0)
        {
          switch (typeOfStat)
            {
            case 0: return TSTEP_AVG;
            case 1: return TSTEP_ACCUM;
            case 2: return TSTEP_MAX;
            case 3: return TSTEP_MIN;
            case 4: return TSTEP_DIFF;
            case 5: return TSTEP_RMS;
            case 6: return TSTEP_SD;
            case 7: return TSTEP_COV;
            case 9: return TSTEP_RATIO;
            case 11: return TSTEP_SUM;
            }
        }

#ifdef HIRLAM_EXTENSIONS
      {
        // Normaly cdo looks in grib for attribute called "stepType", see above.
        // BUT NWP models such as Hirlam and Harmonie 37h1.2, use "timeRangeIndicator" instead!
        // Where for example:       0: for instanteneous fields; 4: for accumulated fields
        //  0:   Forecast product valid at reference time + P1
        //  2:   Product with a valid time ranging between reference time + P1 and reference time + P2
        //  4:   Accumulation (reference time + P1 to reference time + P2)
        //  5:   Difference(reference time + P2 minus reference time + P1) product considered valid at reference time + P2
        // More details on WMO standards:
        //               http://www.wmo.int/pages/prog/www/WDM/Guides/Guide-binary-2.html
        // tsteptype = TSTEP_INSTANT;  // default value for any case
        long timeRangeIND = 0;  // typically 0: for instanteneous fields; 4: for accumulated fields
        int rc = grib_get_long(gh, "timeRangeIndicator", &timeRangeIND);
        if (rc != 0)
          {
            // if ( lprint )
            Warning("Could not get 'stepType' either 'timeRangeIndicator'. Using default!");
            return tsteptype;
          }
        extern int cdiGribUseTimeRangeIndicator;
        cdiGribUseTimeRangeIndicator = 1;
        switch (timeRangeIND)
          {
          case 0: tsteptype = TSTEP_INSTANT; break;
          case 2:
            tsteptype = TSTEP_INSTANT2;
            strcpy(stepType, "instant2");
            break;  // was incorrectly set before into accum
          case 4: tsteptype = TSTEP_ACCUM; break;
          case 5: tsteptype = TSTEP_DIFF; break;
          default:
            if (lprint)
              {
                if (CDI_Debug)
                  Warning("timeRangeIND = %d;  stepType= %s; tsteptype=%d unsupported timeRangeIND at the moment, set to instant!",
                          timeRangeIND, stepType, tsteptype);
                lprint = false;
              }
            break;
          }
        if (CDI_Debug) Warning("timeRangeIND = %d;  stepType= %s; tsteptype=%d", timeRangeIND, stepType, tsteptype);
      }
#endif  // HIRLAM_EXTENSIONS
    }

tsteptypeFound:
  return tsteptype;
}

int
gribGetDatatype(grib_handle *gribHandle)
{
  int datatype;
  if (gribEditionNumber(gribHandle) > 1 && gribCheckString(gribHandle, "packingType", "grid_ieee"))
    {
      datatype = gribCheckLong(gribHandle, "precision", 1) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
    }
  else
    {
      long bitsPerValue;
      datatype = (!grib_get_long(gribHandle, "bitsPerValue", &bitsPerValue) && bitsPerValue > 0 && bitsPerValue <= 32)
                     ? (int) bitsPerValue
                     : CDI_DATATYPE_PACK;
    }
  return datatype;
}

int
gribapiGetParam(grib_handle *gh)
{
  long pdis, pcat, pnum;
  if (gribEditionNumber(gh) <= 1)
    {
      pdis = 255;
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "table2Version", &pcat);
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "indicatorOfParameter", &pnum);
    }
  else
    {
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "discipline", &pdis);
      if (grib_get_long(gh, "parameterCategory", &pcat)) pcat = 0;
      if (grib_get_long(gh, "parameterNumber", &pnum)) pnum = 0;
    }
  return cdiEncodeParam((int) pnum, (int) pcat, (int) pdis);
}

static bool
has_ni(grib_handle *gh)
{
  return (gribGetLong(gh, "Ni") != (long) GRIB_MISSING_LONG);
}

int
gribapiGetGridType(grib_handle *gh)
{
  long gridDefinitionTemplateNumber = gribGetLongDefault(gh, "gridDefinitionTemplateNumber", -1);
  switch (gridDefinitionTemplateNumber)
    {
    case GRIB2_GTYPE_LATLON: return has_ni(gh) ? GRID_LONLAT : GRID_GENERIC;
    case GRIB2_GTYPE_GAUSSIAN: return has_ni(gh) ? GRID_GAUSSIAN : GRID_GAUSSIAN_REDUCED;
    case GRIB2_GTYPE_LATLON_ROT: return GRID_PROJECTION;
    case GRIB2_GTYPE_LCC: return CDI_PROJ_LCC;
    case GRIB2_GTYPE_STERE: return CDI_PROJ_STERE;
    case GRIB2_GTYPE_SPECTRAL: return GRID_SPECTRAL;
    case GRIB2_GTYPE_GME: return GRID_GME;
    case GRIB2_GTYPE_UNSTRUCTURED: return GRID_UNSTRUCTURED;
    default:
      {
        static bool lwarn = true;
        if (lwarn)
          {
            lwarn = false;
            char mesg[256];
            size_t len = sizeof(mesg);
            if (grib_get_string(gh, "gridType", mesg, &len) != 0) mesg[0] = 0;
            Warning("gridDefinitionTemplateNumber %d unsupported (gridType=%s)!", gridDefinitionTemplateNumber, mesg);
          }
      }
    }

  return GRID_GENERIC;
}

static int
gribapiGetIsRotated(grib_handle *gh)
{
  return gribGetLongDefault(gh, "gridDefinitionTemplateNumber", -1) == GRIB2_GTYPE_LATLON_ROT;
}

size_t
gribapiGetGridsize(grib_handle *gh)
{
  size_t gridsize;
  FAIL_ON_GRIB_ERROR(grib_get_size, gh, "values", &gridsize);
  return gridsize;
}

static void
gribapiGetGridGaussianReduced(grib_handle *gh, grid_t *grid, int editionNumber, size_t numberOfPoints)
{
  long lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar);
  grid->np = (int) lpar;

  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nj", &lpar);
  size_t nlat = (size_t) lpar;

  grid->size = numberOfPoints;

  grid->reducedPointsSize = (int) nlat;
  grid->reducedPoints = (int *) Malloc(nlat * sizeof(int));
  long *pl = (long *) Malloc(nlat * sizeof(long));
  size_t dummy = nlat;
  FAIL_ON_GRIB_ERROR(grib_get_long_array, gh, "pl", pl, &dummy);
  for (size_t i = 0; i < nlat; ++i) grid->reducedPoints[i] = (int) pl[i];
  Free(pl);

  grid->y.size = nlat;
  grid->x.inc = 0;
  grid->y.inc = 0;
  grid->x.flag = 0;
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfFirstGridPointInDegrees", &grid->x.first);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfLastGridPointInDegrees", &grid->x.last);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfFirstGridPointInDegrees", &grid->y.first);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfLastGridPointInDegrees", &grid->y.last);

  // FAIL_ON_GRIB_ERROR(grib_get_double, gh, "iDirectionIncrementInDegrees", &grid->x.inc);
  // if ( IS_EQUAL(grid->x.inc, GRIB_MISSING_DOUBLE) ) grid->x.inc = 0;

  if (grid->x.last < grid->x.first)
    {
      if (grid->x.first >= 180.0)
        grid->x.first -= 360.0;
      else
        grid->x.last += 360.0;
    }

  grid->x.flag = 2;

  grid->y.flag = 0;
  // if (IS_NOT_EQUAL(grid->y.first, 0) || IS_NOT_EQUAL(grid->y.last, 0))
  {
    if (grid->y.size > 1)
      {
        if (editionNumber <= 1)
          {
          }
      }
    grid->y.flag = 2;
  }
}

static void
gribapiGetGridRegular(grib_handle *gh, grid_t *grid, int editionNumber, int gridtype, size_t numberOfPoints)
{
  long lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Ni", &lpar);
  size_t nlon = (size_t) lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nj", &lpar);
  size_t nlat = (size_t) lpar;

  if (gridtype == GRID_GAUSSIAN)
    {
      FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfParallelsBetweenAPoleAndTheEquator", &lpar);
      grid->np = (int) lpar;
    }

  if (numberOfPoints != nlon * nlat) Error("numberOfPoints (%zu) and gridSize (%zu) differ!", numberOfPoints, nlon * nlat);

  grid->size = numberOfPoints;
  grid->x.size = nlon;
  grid->y.size = nlat;
  grid->x.inc = 0;
  grid->y.inc = 0;
  grid->x.flag = 0;
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfFirstGridPointInDegrees", &grid->x.first);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "longitudeOfLastGridPointInDegrees", &grid->x.last);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfFirstGridPointInDegrees", &grid->y.first);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "latitudeOfLastGridPointInDegrees", &grid->y.last);
  if (nlon > 1) FAIL_ON_GRIB_ERROR(grib_get_double, gh, "iDirectionIncrementInDegrees", &grid->x.inc);
  if (gridtype == GRID_LONLAT && nlat > 1) FAIL_ON_GRIB_ERROR(grib_get_double, gh, "jDirectionIncrementInDegrees", &grid->y.inc);

  long iscan = 0, jscan = 0;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "iScansNegatively", &iscan);
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "jScansPositively", &jscan);
  if (iscan) grid->x.inc = -grid->x.inc;
  if (!jscan) grid->y.inc = -grid->y.inc;

  if (grid->x.inc < -999 || grid->x.inc > 999) grid->x.inc = 0;
  if (grid->y.inc < -999 || grid->y.inc > 999) grid->y.inc = 0;

  // if ( IS_NOT_EQUAL(grid->x.first, 0) || IS_NOT_EQUAL(grid->x.last, 0) )
  {
    if (grid->x.size > 1)
      {
        // if ( editionNumber <= 1 )
        {
          if (grid->x.last < grid->x.first)
            {
              if (grid->x.first >= 180)
                grid->x.first -= 360;
              else
                grid->x.last += 360;
            }

          // correct xinc if necessary
          if (IS_EQUAL(grid->x.first, 0) && grid->x.last > 354 && grid->x.last < 360)
            {
              double xinc = 360. / grid->x.size;
              if (fabs(grid->x.inc - xinc) > 0.0)
                {
                  grid->x.inc = xinc;
                  if (CDI_Debug) Message("set xinc to %g", grid->x.inc);
                }
            }
        }
      }
    grid->x.flag = 2;
  }

  grid->y.flag = 0;
  // if ( IS_NOT_EQUAL(grid->y.first, 0) || IS_NOT_EQUAL(grid->y.last, 0) )
  {
    if (grid->y.size > 1)
      {
        if (editionNumber <= 1)
          {
          }
      }
    grid->y.flag = 2;
  }
}

static void
gribapiGetGridProj(grib_handle *gh, grid_t *grid, size_t numberOfPoints)
{
  long lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Nx", &lpar);
  size_t nlon = (size_t) lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "Ny", &lpar);
  size_t nlat = (size_t) lpar;

  if (numberOfPoints != nlon * nlat) Error("numberOfPoints (%zu) and gridSize (%zu) differ!", numberOfPoints, nlon * nlat);

  grid->size = numberOfPoints;
  grid->x.size = nlon;
  grid->y.size = nlat;

  double xinc, yinc;
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "DxInMetres", &xinc);
  FAIL_ON_GRIB_ERROR(grib_get_double, gh, "DyInMetres", &yinc);

  grid->x.first = 0;
  grid->x.last = 0;
  grid->x.inc = xinc;
  grid->y.first = 0;
  grid->y.last = 0;
  grid->y.inc = yinc;
  grid->x.flag = 2;
  grid->y.flag = 2;
}

static void
gribapiGetGridSpectral(grib_handle *gh, grid_t *grid, size_t datasize)
{
  size_t len = 256;
  char typeOfPacking[256];
  FAIL_ON_GRIB_ERROR(grib_get_string, gh, "packingType", typeOfPacking, &len);
  grid->lcomplex = 0;
  if (strncmp(typeOfPacking, "spectral_complex", len) == 0) grid->lcomplex = 1;

  grid->size = datasize;

  long lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "J", &lpar);
  grid->trunc = (int) lpar;
}

static void
gribapiGetGridGME(grib_handle *gh, grid_t *grid, size_t numberOfPoints)
{
  grid->size = numberOfPoints;

  long lpar;
  if (grib_get_long(gh, "nd", &lpar) == 0) grid->gme.nd = (int) lpar;
  if (grib_get_long(gh, "Ni", &lpar) == 0) grid->gme.ni = (int) lpar;
  if (grib_get_long(gh, "n2", &lpar) == 0) grid->gme.ni2 = (int) lpar;
  if (grib_get_long(gh, "n3", &lpar) == 0) grid->gme.ni3 = (int) lpar;
}

static void
gribapiGetGridUnstructured(grib_handle *gh, grid_t *grid, size_t numberOfPoints)
{
  unsigned char uuid[CDI_UUID_SIZE];
  /*
    char reference_link[8192];
    size_t len = sizeof(reference_link);
    reference_link[0] = 0;
  */
  grid->size = numberOfPoints;

  long lpar;
  if (grib_get_long(gh, "numberOfGridUsed", &lpar) == 0)
    {
      cdiDefVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDUSED, (int) lpar);
      if (grib_get_long(gh, "numberOfGridInReference", &lpar) == 0)
        cdiDefVarKeyInt(&grid->keys, CDI_KEY_NUMBEROFGRIDINREFERENCE, (int) lpar);
      /*
        if ( grib_get_string(gh, "gridDescriptionFile", reference_link, &len) == 0 )
        {
        if ( strncmp(reference_link, "file://", 7) == 0 )
        grid->reference = strdupx(reference_link);
        }
      */
      size_t len = (size_t) CDI_UUID_SIZE;
      if (grib_get_bytes(gh, "uuidOfHGrid", uuid, &len) == 0) cdiDefVarKeyBytes(&grid->keys, CDI_KEY_UUID, uuid, CDI_UUID_SIZE);
    }
}

static void
gribapiGetGridGeneric(grib_handle *gh, grid_t *grid, size_t numberOfPoints)
{
  long lpar;
  size_t nlon = (grib_get_long(gh, "Ni", &lpar) == 0) ? (size_t) lpar : 0;
  size_t nlat = (grib_get_long(gh, "Nj", &lpar) == 0) ? (size_t) lpar : 0;

  grid->size = numberOfPoints;

  bool lgeneric = (nlon > 0 && nlat > 0 && nlon * nlat == numberOfPoints);
  grid->x.size = lgeneric ? nlon : 0;
  grid->y.size = lgeneric ? nlat : 0;
}

// TODO: Simplify by use of the convenience functions (gribGetLong(), gribGetLongDefault(), etc.).
bool
gribapiGetGrid(grib_handle *gh, grid_t *grid)
{
  bool uvRelativeToGrid = false;
  long editionNumber = gribEditionNumber(gh);
  int gridtype = gribapiGetGridType(gh);
  int projtype = (gridtype == GRID_PROJECTION && gribapiGetIsRotated(gh)) ? CDI_PROJ_RLL : CDI_UNDEFID;
  if (gridtype == CDI_PROJ_LCC || gridtype == CDI_PROJ_STERE)
    {
      projtype = gridtype;
      gridtype = GRID_PROJECTION;
    }
  /*
  if ( streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED )
    {
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = 2*ISEC2_NumLat;
      ISEC4_NumValues = ISEC2_NumLon*ISEC2_NumLat;
    }
  */
  grid_init(grid);
  cdiGridTypeInit(grid, gridtype, 0);

  size_t datasize;
  FAIL_ON_GRIB_ERROR(grib_get_size, gh, "values", &datasize);
  long lpar;
  FAIL_ON_GRIB_ERROR(grib_get_long, gh, "numberOfPoints", &lpar);
  size_t numberOfPoints = (size_t) lpar;

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_RLL)
    {
      gribapiGetGridRegular(gh, grid, editionNumber, gridtype, numberOfPoints);
    }
  else if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      gribapiGetGridGaussianReduced(gh, grid, editionNumber, numberOfPoints);
    }
  else if (projtype == CDI_PROJ_LCC)
    {
      gribapiGetGridProj(gh, grid, numberOfPoints);
    }
  else if (projtype == CDI_PROJ_STERE)
    {
      gribapiGetGridProj(gh, grid, numberOfPoints);
    }
  else if (gridtype == GRID_SPECTRAL)
    {
      gribapiGetGridSpectral(gh, grid, datasize);
    }
  else if (gridtype == GRID_GME)
    {
      gribapiGetGridGME(gh, grid, numberOfPoints);
    }
  else if (gridtype == GRID_UNSTRUCTURED)
    {
      gribapiGetGridUnstructured(gh, grid, numberOfPoints);
    }
  else if (gridtype == GRID_GENERIC)
    {
      gribapiGetGridGeneric(gh, grid, numberOfPoints);
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || projtype == CDI_PROJ_RLL || projtype == CDI_PROJ_LCC)
    {
      long temp = 0;
      GRIB_CHECK(grib_get_long(gh, "uvRelativeToGrid", &temp), 0);
      assert(temp == 0 || temp == 1);
      uvRelativeToGrid = (bool) temp;
    }

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION)
    {
      long iScansNegatively, jScansPositively, jPointsAreConsecutive;
      GRIB_CHECK(grib_get_long(gh, "iScansNegatively", &iScansNegatively), 0);
      GRIB_CHECK(grib_get_long(gh, "jScansPositively", &jScansPositively), 0);
      GRIB_CHECK(grib_get_long(gh, "jPointsAreConsecutive", &jPointsAreConsecutive), 0);

      int scanningMode = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
      cdiDefVarKeyInt(&grid->keys, CDI_KEY_SCANNINGMODE, scanningMode);
      /* scanningMode  = 128 * iScansNegatively + 64 * jScansPositively + 32 * jPointsAreConsecutive;
                   64  = 128 * 0                + 64 *        1         + 32 * 0
                   00  = 128 * 0                + 64 *        0         + 32 * 0
                   96  = 128 * 0                + 64 *        1         + 32 * 1
         Default / implicit scanning mode is 64:
                            i and j scan positively, i points are consecutive (row-major)        */
#ifdef HIRLAM_EXTENSIONS
      if (cdiDebugExt >= 30 && gribEditionNumber(gh) <= 1)
        {
          //  indicatorOfParameter=33,indicatorOfTypeOfLevel=105,level
          long paramId, levelTypeId, levelId;
          GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &paramId), 0);
          GRIB_CHECK(grib_get_long(gh, "indicatorOfTypeOfLevel", &levelTypeId), 0);
          GRIB_CHECK(grib_get_long(gh, "level", &levelId), 0);
          Message("(param,ltype,level) = (%3d,%3d,%4d); Scanning mode = %02d -> bits:(%1d.%1d.%1d)*32;  uvRelativeToGrid = %02d",
                  (int) paramId, (int) levelTypeId, (int) levelId, scanningMode, jPointsAreConsecutive, jScansPositively,
                  iScansNegatively, uvRelativeToGrid);
        }
#endif  // HIRLAM_EXTENSIONS
    }

  grid->type = gridtype;
  grid->projtype = projtype;

  return uvRelativeToGrid;
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
