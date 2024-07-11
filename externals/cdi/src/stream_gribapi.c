#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBGRIB_API

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "gribapi_utilities.h"
#include "stream_scan.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "varscan.h"
#include "vlist.h"
#include "subtype.h"
#include "cdi_uuid.h"

#include "cgribex.h" /* gribGetSize, gribRead, gribGetZip, GRIB1_LTYPE_99 */
#include "gribapi.h"

#include <grib_api.h>

extern int CDI_Inventory_Mode;

static const var_tile_t dummy_tiles = { 0, -1, -1, -1, -1, -1 };

typedef struct
{
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
  size_t gridsize;
  char name[32];
  VarScanKeys scanKeys;
  var_tile_t tiles;
} compvar2_t;

static int
gribapi_get_zaxis_type(long editionNumber, int grib_ltype)
{
  return (editionNumber <= 1) ? grib1ltypeToZaxisType(grib_ltype) : grib2ltypeToZaxisType(grib_ltype);
}

static int
get_timeunits(long unitsOfTime)
{
  switch (unitsOfTime)
    {
    case 13: return TUNIT_SECOND;
    case 0: return TUNIT_MINUTE;
    case 1: return TUNIT_HOUR;
    case 10: return TUNIT_3HOURS;
    case 11: return TUNIT_6HOURS;
    case 12: return TUNIT_12HOURS;
    case 2: return TUNIT_DAY;
    }

  return TUNIT_HOUR;
}

static double
timeunits_factor(int timeUnits1, int timeUnits2)
{
  if (timeUnits2 == TUNIT_HOUR)
    {
      switch (timeUnits1)
        {
        case TUNIT_SECOND: return 3600;
        case TUNIT_MINUTE: return 60;
        case TUNIT_HOUR: return 1;
        case TUNIT_3HOURS: return 1.0 / 3;
        case TUNIT_6HOURS: return 1.0 / 6;
        case TUNIT_12HOURS: return 1.0 / 12;
        case TUNIT_DAY: return 1.0 / 24;
        }
    }

  return 1;
}

static int
gribapi_get_timeunits(grib_handle *gh)
{
  long unitsOfTime = -1;
  grib_get_long(gh, "indicatorOfUnitOfTimeRange", &unitsOfTime);

  GRIB_CHECK(my_grib_set_long(gh, "stepUnits", unitsOfTime), 0);

  return get_timeunits(unitsOfTime);
}

static void
gribapiGetSteps(grib_handle *gh, int timeunits, int *startStep, int *endStep)
{
  long unitsOfTime;
  int status = grib_get_long(gh, "stepUnits", &unitsOfTime);
  int timeunits2 = (status == 0) ? get_timeunits(unitsOfTime) : timeunits;
  // timeunits2 = gribapi_get_timeunits(gh);

  long lpar;
  status = grib_get_long(gh, "forecastTime", &lpar);
  if (status == 0)
    *startStep = (int) lpar;
  else
    {
      status = grib_get_long(gh, "startStep", &lpar);
      if (status == 0) *startStep = (int) (((double) lpar * timeunits_factor(timeunits, timeunits2)) + 0.5);
    }

  *endStep = *startStep;
  status = grib_get_long(gh, "endStep", &lpar);
  if (status == 0) *endStep = (int) (((double) lpar * timeunits_factor(timeunits, timeunits2)) + 0.5);
  // printf("%d %d %d %d %d %g\n", *startStep, *endStep, lpar, timeunits, timeunits2, timeunits_factor(timeunits, timeunits2));
}

static CdiDateTime
gribapiGetDataDateTime(grib_handle *gh)
{
  long date;
  GRIB_CHECK(grib_get_long(gh, "dataDate", &date), 0);

  long hour, minute, second;
  GRIB_CHECK(grib_get_long(gh, "hour", &hour), 0);
  GRIB_CHECK(grib_get_long(gh, "minute", &minute), 0);
  GRIB_CHECK(grib_get_long(gh, "second", &second), 0);

  CdiDateTime dt;
  dt.date = cdiDate_set(date);
  dt.time.hour = hour;
  dt.time.minute = minute;
  dt.time.second = second;
  dt.time.ms = 0;

  return dt;
}

static void
gribapiSetDataDateTime(grib_handle *gh, CdiDateTime dataDateTime)
{
  long dataDate = (long) cdiDate_get(dataDateTime.date);
  GRIB_CHECK(my_grib_set_long(gh, "dataDate", (long) dataDate), 0);

  int hour, minute, second, ms;
  cdiTime_decode(dataDateTime.time, &hour, &minute, &second, &ms);
  GRIB_CHECK(my_grib_set_long(gh, "hour", hour), 0);
  GRIB_CHECK(my_grib_set_long(gh, "minute", minute), 0);
  GRIB_CHECK(my_grib_set_long(gh, "second", second), 0);
}

static int
gribapi_get_timeunits_factor(int timeUnits)
{
  static bool lprint = true;
  switch (timeUnits)
    {
    case TUNIT_SECOND: return 1;
    case TUNIT_MINUTE: return 60;
    case TUNIT_HOUR: return 3600;
    case TUNIT_3HOURS: return 10800;
    case TUNIT_6HOURS: return 21600;
    case TUNIT_12HOURS: return 43200;
    case TUNIT_DAY: return 86400;
    default:
      if (lprint)
        {
          Warning("Time unit %d unsupported", timeUnits);
          lprint = false;
        }
    }

  return 0;
}

static CdiDateTime
gribapiGetValidityDateTime(grib_handle *gh, CdiDateTime *sDateTime)
{
  CdiDateTime vDateTime;
  cdiDateTime_init(sDateTime);

  long sigofrtime = 3;
  if (gribEditionNumber(gh) > 1)
    GRIB_CHECK(grib_get_long(gh, "significanceOfReferenceTime", &sigofrtime), 0);
  else
    GRIB_CHECK(grib_get_long(gh, "timeRangeIndicator", &sigofrtime), 0);

  if (sigofrtime
      == 3)  // XXX: This looks like a bug to me, because timeRangeIndicator == 3 does not seem to have the same meaning as
             // significanceOfReferenceTime == 3. I would recommend replacing this condition with `if(!gribapiTimeIsFC())`.
    {
      vDateTime = gribapiGetDataDateTime(gh);
    }
  else
    {
      CdiDateTime rDateTime = gribapiGetDataDateTime(gh);

      int timeUnits = gribapi_get_timeunits(gh);
      int startStep = 0, endStep = 0;
      gribapiGetSteps(gh, timeUnits, &startStep, &endStep);

      if (rDateTime.date.day > 0)
        {
          extern int CGRIBEX_grib_calendar;
          JulianDate julianDate = julianDate_encode(CGRIBEX_grib_calendar, rDateTime);

          int64_t timeUnitFactor = gribapi_get_timeunits_factor(timeUnits);

          // if (startStep > 0)
          {
            JulianDate julianDate2 = julianDate_add_seconds(julianDate, timeUnitFactor * startStep);
            *sDateTime = julianDate_decode(CGRIBEX_grib_calendar, julianDate2);
          }

          rDateTime = julianDate_decode(CGRIBEX_grib_calendar, julianDate_add_seconds(julianDate, timeUnitFactor * endStep));
        }

      vDateTime = rDateTime;
    }

  return vDateTime;
}

static void
grib1_get_level(grib_handle *gh, int *leveltype, int *lbounds, int *level1, int *level2)
{
  *leveltype = 0;
  *lbounds = 0;
  *level1 = 0;
  *level2 = 0;

  long lpar;
  if (!grib_get_long(gh, "indicatorOfTypeOfLevel", &lpar))  // 1 byte
    {
      *leveltype = (int) lpar;

      switch (*leveltype)
        {
        case GRIB1_LTYPE_SIGMA_LAYER:
        case GRIB1_LTYPE_HYBRID_LAYER:
        case GRIB1_LTYPE_LANDDEPTH_LAYER:
          {
            *lbounds = 1;
            break;
          }
        }

      if (*lbounds)
        {
          GRIB_CHECK(grib_get_long(gh, "topLevel", &lpar), 0);  // 1 byte
          if (lpar == GRIB_MISSING_LONG) lpar = 0;
          *level1 = (int) lpar;
          GRIB_CHECK(grib_get_long(gh, "bottomLevel", &lpar), 0);  // 1 byte
          if (lpar == GRIB_MISSING_LONG) lpar = 0;
          *level2 = (int) lpar;
        }
      else
        {
          double dlevel;
          GRIB_CHECK(grib_get_double(gh, "level", &dlevel), 0);  // 2 byte
          if (*leveltype == GRIB1_LTYPE_ISOBARIC) dlevel *= 100;
          if (dlevel < -2.e9 || dlevel > 2.e9) dlevel = 0;
          if (*leveltype == GRIB1_LTYPE_99 || *leveltype == GRIB1_LTYPE_ISOBARIC_PA) *leveltype = GRIB1_LTYPE_ISOBARIC;

          *level1 = (int) dlevel;
          *level2 = 0;
        }
    }
}

static double
grib2ScaleFactor(long factor)
{
  switch (factor)
    {
    case GRIB_MISSING_LONG: return 1;
    case -9: return 1000000000;
    case -8: return 100000000;
    case -7: return 10000000;
    case -6: return 1000000;
    case -5: return 100000;
    case -4: return 10000;
    case -3: return 1000;
    case -2: return 100;
    case -1: return 10;
    case 0: return 1;
    case 1: return 0.1;
    case 2: return 0.01;
    case 3: return 0.001;
    case 4: return 0.0001;
    case 5: return 0.00001;
    case 6: return 0.000001;
    case 7: return 0.0000001;
    case 8: return 0.00000001;
    case 9: return 0.000000001;
    default: return 0;
    }
}

static int
calc_level(int level_sf, long factor, long level)
{
  double result = 0;
  if (level != GRIB_MISSING_LONG) result = (double) level * grib2ScaleFactor(factor);
  if (level_sf) result *= level_sf;
  return (int) result;
}

static void
grib2_get_level(grib_handle *gh, int *leveltype1, int *leveltype2, int *lbounds, int *level1, int *level2, int *level_sf,
                int *level_unit)
{
  *leveltype1 = 0;
  *leveltype2 = -1;
  *lbounds = 0;
  *level1 = 0;
  *level2 = 0;
  *level_sf = 0;
  *level_unit = 0;

  long lpar;
  int status = grib_get_long(gh, "typeOfFirstFixedSurface", &lpar);  // 1 byte
  if (status == 0)
    {
      *leveltype1 = (int) lpar;

      status = grib_get_long(gh, "typeOfSecondFixedSurface", &lpar);  // 1 byte
      /* FIXME: assert(lpar >= INT_MIN && lpar <= INT_MAX) */
      if (status == 0) *leveltype2 = (int) lpar;

      if (*leveltype1 != 255 && *leveltype2 != 255 && *leveltype2 > 0) *lbounds = 1;
      switch (*leveltype1)
        {
        case GRIB2_LTYPE_REFERENCE:
          if (*leveltype2 == 1) *lbounds = 0;
          break;
        case GRIB2_LTYPE_LANDDEPTH:
          *level_sf = 1000;
          *level_unit = CDI_UNIT_M;
          break;
        case GRIB2_LTYPE_ISOBARIC:
          *level_sf = 1000;
          *level_unit = CDI_UNIT_PA;
          break;
        case GRIB2_LTYPE_SIGMA:
          *level_sf = 1000;
          *level_unit = 0;
          break;
        }

      long factor, llevel;
      GRIB_CHECK(grib_get_long(gh, "scaleFactorOfFirstFixedSurface", &factor), 0);  // 1 byte
      GRIB_CHECK(grib_get_long(gh, "scaledValueOfFirstFixedSurface", &llevel), 0);  // 4 byte
      *level1 = calc_level(*level_sf, factor, llevel);

      if (*lbounds)
        {
          GRIB_CHECK(grib_get_long(gh, "scaleFactorOfSecondFixedSurface", &factor), 0);  // 1 byte
          GRIB_CHECK(grib_get_long(gh, "scaledValueOfSecondFixedSurface", &llevel), 0);  // 4 byte
          *level2 = calc_level(*level_sf, factor, llevel);
        }
    }
}

static void
grib_get_level(grib_handle *gh, int *leveltype1, int *leveltype2, int *lbounds, int *level1, int *level2, int *level_sf,
               int *level_unit, var_tile_t *tiles)
{
  if (gribEditionNumber(gh) <= 1)
    {
      grib1_get_level(gh, leveltype1, lbounds, level1, level2);
      *leveltype2 = -1;
      *level_sf = 0;
      *level_unit = 0;
    }
  else
    {
      grib2_get_level(gh, leveltype1, leveltype2, lbounds, level1, level2, level_sf, level_unit);

      // read in tiles attributes (if there are any)
      tiles->tileindex = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEINDEX], 0);
      tiles->totalno_of_tileattr_pairs
          = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS], -1);
      tiles->tileClassification = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILE_CLASSIFICATION], -1);
      tiles->numberOfTiles = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_TILES], -1);
      tiles->numberOfAttributes = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_NUMBER_OF_ATTR], -1);
      tiles->attribute = (int) gribGetLongDefault(gh, cdiSubtypeAttributeName[SUBTYPE_ATT_TILEATTRIBUTE], -1);
    }
}

static void
gribapi_get_string(grib_handle *gh, const char *key, char *string, size_t length)
{
  string[0] = 0;

  int ret = grib_get_string(gh, key, string, &length);
  if (ret != 0)
    {
      fprintf(stderr, "grib_get_string(gh, \"%s\", ...) failed!\n", key);
      GRIB_CHECK(ret, 0);
    }
  // clang-format off
  if      (length == 8 && memcmp(string, "unknown", length) == 0) string[0] = 0;
  else if (length == 2 && memcmp(string, "~", length) == 0)       string[0] = 0;
  // clang-format on
}

static void
param_to_name(int param, char *name)
{
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  if (pdis == 255)
    {
      snprintf(name, 256, "code%d", pnum);
    }
  else
    {
      snprintf(name, 256, "param%d.%d.%d", pnum, pcat, pdis);
    }
}

static int
gribapiGetEnsembleInfo(grib_handle *gh, long *numberOfForecastsInEnsemble, long *perturbationNumber, long *typeOfEnsembleForecast)
{
  int status = 0;
  if (grib_get_long(gh, "numberOfForecastsInEnsemble", numberOfForecastsInEnsemble) == 0)
    {
      if (*numberOfForecastsInEnsemble > 0) status = 1;
      GRIB_CHECK(grib_get_long(gh, "perturbationNumber", perturbationNumber), 0);
      grib_get_long(gh, "typeOfEnsembleForecast", typeOfEnsembleForecast);
    }

  if (status == 0)
    {
      *numberOfForecastsInEnsemble = -1;
      *perturbationNumber = -1;
      *typeOfEnsembleForecast = -1;
    }

  return status;
}

static VarScanKeys
gribapiGetScanKeys(grib_handle *gh)
{
  VarScanKeys scanKeys;
  varScanKeysInit(&scanKeys);

  long numberOfForecastsInEnsemble = -1, perturbationNumber = -1, typeOfEnsembleForecast = -1;
  gribapiGetEnsembleInfo(gh, &numberOfForecastsInEnsemble, &perturbationNumber, &typeOfEnsembleForecast);
  scanKeys.perturbationNumber = (short) perturbationNumber;

  long typeOfGeneratingProcess = 0;
  if (grib_get_long(gh, "typeOfGeneratingProcess", &typeOfGeneratingProcess) == 0)
    scanKeys.typeOfGeneratingProcess = (short) typeOfGeneratingProcess;

  return scanKeys;
}

static void
gribapiGetNameKeys(grib_handle *gh, int varID)
{
  char string[CDI_MAX_NAME];

  size_t vlen = CDI_MAX_NAME;
  gribapi_get_string(gh, "name", string, vlen);  // longname
  if (string[0]) varDefKeyString(varID, CDI_KEY_LONGNAME, string);

  gribapi_get_string(gh, "units", string, vlen);
  if (string[0]) varDefKeyString(varID, CDI_KEY_UNITS, string);

  string[0] = 0;
  int status = grib_get_string(gh, "cfName", string, &vlen);
  if (status != 0 || vlen <= 1 || strncmp(string, "unknown", 7) == 0) string[0] = 0;
  if (string[0]) varDefKeyString(varID, CDI_KEY_STDNAME, string);
}

static void
gribapiGetKeys(grib_handle *gh, int varID)
{
  long tablesVersion = 0;
  if (grib_get_long(gh, "tablesVersion", &tablesVersion) == 0) varDefKeyInt(varID, CDI_KEY_TABLESVERSION, (int) tablesVersion);

  long localTablesVersion = 0;
  if (grib_get_long(gh, "localTablesVersion", &localTablesVersion) == 0)
    varDefKeyInt(varID, CDI_KEY_LOCALTABLESVERSION, (int) localTablesVersion);

  long typeOfGeneratingProcess = 0;
  if (grib_get_long(gh, "typeOfGeneratingProcess", &typeOfGeneratingProcess) == 0)
    varDefKeyInt(varID, CDI_KEY_TYPEOFGENERATINGPROCESS, (int) typeOfGeneratingProcess);

  long productDefinitionTemplate = 0;
  if (grib_get_long(gh, "productDefinitionTemplateNumber", &productDefinitionTemplate) == 0)
    varDefKeyInt(varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, (int) productDefinitionTemplate);

  long typeOfProcessedData = 0;
  if (grib_get_long(gh, "typeOfProcessedData", &typeOfProcessedData) == 0)
    varDefKeyInt(varID, CDI_KEY_TYPEOFPROCESSEDDATA, (int) typeOfProcessedData);

  long shapeOfTheEarth = 0;
  if (grib_get_long(gh, "shapeOfTheEarth", &shapeOfTheEarth) == 0)
    varDefKeyInt(varID, CDI_KEY_SHAPEOFTHEEARTH, (int) shapeOfTheEarth);

  long backgroundProcess = 0;
  if (grib_get_long(gh, "backgroundProcess", &backgroundProcess) == 0)
    varDefKeyInt(varID, CDI_KEY_BACKGROUNDPROCESS, (int) backgroundProcess);

  long typeOfTimeIncrement = 0;
  if (grib_get_long(gh, "typeOfTimeIncrement", &typeOfTimeIncrement) == 0)
    varDefKeyInt(varID, CDI_KEY_TYPEOFTIMEINCREMENT, (int) typeOfTimeIncrement);
  /*
  long constituentType = 0;
  if ( grib_get_long(gh, "constituentType", &constituentType) == 0 )
    varDefKeyInt(varID, CDI_KEY_CONSTITUENTTYPE, (int) constituentType);
  */

  /*
    Get the ensemble information from the grib-2 Tables and update the intermediate datastructure.
    Further update to the "vlist" is handled in the same way as for GRIB-1 by "cdi_generate_vars"
  */
  long numberOfForecastsInEnsemble = -1, perturbationNumber = -1, typeOfEnsembleForecast = -1;
  gribapiGetEnsembleInfo(gh, &numberOfForecastsInEnsemble, &perturbationNumber, &typeOfEnsembleForecast);
  if (numberOfForecastsInEnsemble > 0)
    {
      varDefKeyInt(varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, (int) numberOfForecastsInEnsemble);
      varDefKeyInt(varID, CDI_KEY_PERTURBATIONNUMBER, (int) perturbationNumber);
      if (typeOfEnsembleForecast != -1) varDefKeyInt(varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, (int) typeOfEnsembleForecast);
    }

  long section2Length = 0;
  int status = grib_get_long(gh, "section2Length", &section2Length);
  if (status == 0 && section2Length > 0)
    {
      long grib2LocalSectionNumber;
      long mpimType, mpimClass, mpimUser;
      status = grib_get_long(gh, "grib2LocalSectionNumber", &grib2LocalSectionNumber);
      if (status == 0)
        {
          size_t section2PaddingLength = 0;
          status = grib_get_size(gh, "section2Padding", &section2PaddingLength);
          if (status == 0 && section2PaddingLength > 0)
            {
              varDefKeyInt(varID, CDI_KEY_GRIB2LOCALSECTIONNUMBER, (int) grib2LocalSectionNumber);
              varDefKeyInt(varID, CDI_KEY_SECTION2PADDINGLENGTH, (int) section2PaddingLength);
              unsigned char *section2Padding = (unsigned char *) Malloc(section2PaddingLength);
              grib_get_bytes(gh, "section2Padding", section2Padding, &section2PaddingLength);
              varDefKeyBytes(varID, CDI_KEY_SECTION2PADDING, section2Padding, (int) section2PaddingLength);
              Free(section2Padding);
            }
          else if (grib_get_long(gh, "mpimType", &mpimType) == 0 && grib_get_long(gh, "mpimClass", &mpimClass) == 0
                   && grib_get_long(gh, "mpimUser", &mpimUser) == 0)
            {
              varDefKeyInt(varID, CDI_KEY_MPIMTYPE, mpimType);
              varDefKeyInt(varID, CDI_KEY_MPIMCLASS, mpimClass);
              varDefKeyInt(varID, CDI_KEY_MPIMUSER, mpimUser);

              size_t revNumLen = 20;
              unsigned char revNumber[revNumLen];
              if (grib_get_bytes(gh, "revNumber", revNumber, &revNumLen) == 0)
                varDefKeyBytes(varID, CDI_KEY_REVNUMBER, revNumber, (int) revNumLen);

              long revStatus;
              grib_get_long(gh, "revStatus", &revStatus);
              varDefKeyInt(varID, CDI_KEY_REVSTATUS, revStatus);
            }
        }
    }
}

static void
gribapiDefProjRLL(grib_handle *gh, int gridID)
{
  double xpole = 0, ypole = 0, angle = 0;
  grib_get_double(gh, "latitudeOfSouthernPoleInDegrees", &ypole);
  grib_get_double(gh, "longitudeOfSouthernPoleInDegrees", &xpole);
  grib_get_double(gh, "angleOfRotation", &angle);
  xpole -= 180;
  if (fabs(ypole) > 0) ypole = -ypole;  // change from south to north pole
  if (fabs(angle) > 0) angle = -angle;

  gridDefParamRLL(gridID, xpole, ypole, angle);
}

static void
decode_shapeOfTheEarth(grib_handle *gh, struct CDI_GridProjParams *gpp)
{
  long shapeOfTheEarth = 0;
  grib_get_long(gh, "shapeOfTheEarth", &shapeOfTheEarth);

  double radiusOfTheEarth = 6367470.0;
  if (shapeOfTheEarth == 1) grib_get_double(gh, "radiusOfTheEarth", &radiusOfTheEarth);

  // clang-format off
  switch (shapeOfTheEarth)
    {
    case 0:   gpp->a = radiusOfTheEarth; break;
    case 1:   gpp->a = radiusOfTheEarth; break;
    case 2:   gpp->a = 6378160.0; gpp->b = 6356775.0;   gpp->rf = 297.0; break;
    case 4:   gpp->a = 6378137.0; gpp->b = 6356752.314; gpp->rf = 298.257222101; break;
    case 6:   gpp->a = 6371229.0; break;
    case 8:   gpp->a = 6371200.0; break;
    default:  gpp->a = radiusOfTheEarth;
    }
  // clang-format on
}

static void
gribapiDefProjLCC(grib_handle *gh, int gridID)
{
  struct CDI_GridProjParams gpp;
  gridProjParamsInit(&gpp);

  decode_shapeOfTheEarth(gh, &gpp);

  long projflag = 0;
  grib_get_double(gh, "longitudeOfFirstGridPointInDegrees", &gpp.xval_0);
  grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &gpp.yval_0);
  grib_get_double(gh, "longitudeOfSouthernPoleInDegrees", &gpp.x_SP);
  grib_get_double(gh, "latitudeOfSouthernPoleInDegrees", &gpp.y_SP);
  grib_get_double(gh, "LoVInDegrees", &gpp.lon_0);
  grib_get_double(gh, "Latin1InDegrees", &gpp.lat_1);
  grib_get_double(gh, "Latin2InDegrees", &gpp.lat_2);
  grib_get_long(gh, "projectionCentreFlag", &projflag);
  bool isSouthPole = gribbyte_get_bit((int) projflag, 1);
  if (isSouthPole)
    {
      gpp.lat_1 = -gpp.lat_1;
      gpp.lat_2 = -gpp.lat_2;
    }

  gpp.lat_0 = gpp.lat_2;

  if (proj_lonlat_to_lcc_func)
    {
      double x_0 = gpp.xval_0, y_0 = gpp.yval_0;
      proj_lonlat_to_lcc_func(gpp, (size_t) 1, &x_0, &y_0);
      if (IS_NOT_EQUAL(x_0, gpp.mv) && IS_NOT_EQUAL(y_0, gpp.mv))
        {
          gpp.x_0 = -x_0;
          gpp.y_0 = -y_0;
        }
    }

  gridDefParamsLCC(gridID, gpp);
}

static void
gribapiDefProjSTERE(grib_handle *gh, int gridID)
{
  struct CDI_GridProjParams gpp;
  gridProjParamsInit(&gpp);

  decode_shapeOfTheEarth(gh, &gpp);

  grib_get_double(gh, "longitudeOfFirstGridPointInDegrees", &gpp.xval_0);
  grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &gpp.yval_0);
  grib_get_double(gh, "LaDInDegrees", &gpp.lat_1);
  grib_get_double(gh, "orientationOfTheGridInDegrees", &gpp.lon_0);

  long southPoleOnProjectionPlane;
  grib_get_long(gh, "southPoleOnProjectionPlane", &southPoleOnProjectionPlane);
  gpp.lat_0 = southPoleOnProjectionPlane ? -90.0 : 90.0;

  if (proj_lonlat_to_stere_func)
    {
      double x_0 = gpp.xval_0, y_0 = gpp.yval_0;
      proj_lonlat_to_stere_func(gpp, (size_t) 1, &x_0, &y_0);
      if (IS_NOT_EQUAL(x_0, gpp.mv) && IS_NOT_EQUAL(y_0, gpp.mv))
        {
          gpp.x_0 = -x_0;
          gpp.y_0 = -y_0;
        }
    }

  gridDefParamsSTERE(gridID, gpp);
}

static void
gribapiDefProjHEALPIX(grib_handle *gh, int gridID)
{
  struct CDI_GridProjParams gpp;
  gridProjParamsInit(&gpp);

  decode_shapeOfTheEarth(gh, &gpp);

  long lval = -1;
  grib_get_long(gh, "Nside", &lval);
  gpp.nside = (int) lval;
  lval = -1;
  grib_get_long(gh, "ordering", &lval);
  gpp.order = (int) lval;

  gridDefParamsHEALPIX(gridID, gpp);
}

static void
gribapiAddRecord(stream_t *streamptr, int param, grib_handle *gh, size_t recsize, off_t position, int datatype, int comptype,
                 const char *varname, int leveltype1, int leveltype2, int lbounds, int level1, int level2, int level_sf,
                 int level_unit, VarScanKeys *scanKeys, const var_tile_t *tiles, bool lread_additional_keys, int fdbItemIndex)
{
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int recID = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  int tsteptype = gribapiGetTsteptype(gh);

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype1);

  record->size = recsize;
  record->position = position;
  record->param = param;
  record->ilevel = level1;
  record->ilevel2 = level2;
  record->ltype = leveltype1;
  record->tsteptype = (short) tsteptype;
  record->gridsize = gribapiGetGridsize(gh);
  record->scanKeys = *scanKeys;
  record->tiles = tiles ? *tiles : dummy_tiles;
#ifdef HAVE_LIBFDB5
  record->fdbItemIndex = fdbItemIndex;
#else
  (void) fdbItemIndex;
#endif

  strncpy(record->varname, varname, sizeof(record->varname) - 1);
  record->varname[sizeof(record->varname) - 1] = 0;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  CdiQuery *query = streamptr->query;
  if (query && cdiQueryName(query, varname) < 0)
    {
      record->used = false;
      return;
    }

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  grid_t *gridptr = (grid_t *) Malloc(sizeof(*gridptr));
  bool uvRelativeToGrid = gribapiGetGrid(gh, gridptr);

  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, gridptr, 0);
  int gridID = gridAdded.Id;
  // clang-format off
  if (!gridAdded.isNew)
    {
      grid_free(gridptr);
      Free(gridptr);
    }
  else if (gridptr->projtype == CDI_PROJ_RLL)     gribapiDefProjRLL(gh, gridID);
  else if (gridptr->projtype == CDI_PROJ_LCC)     gribapiDefProjLCC(gh, gridID);
  else if (gridptr->projtype == CDI_PROJ_STERE)   gribapiDefProjSTERE(gh, gridID);
  else if (gridptr->projtype == CDI_PROJ_HEALPIX) gribapiDefProjHEALPIX(gh, gridID);
  // clang-format on

  int zaxistype = gribapi_get_zaxis_type(gribEditionNumber(gh), leveltype1);

  switch (zaxistype)
    {
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        long lpar;
        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        /* FIXME: assert(lpar >= 0) */
        size_t vctsize = (size_t) lpar;
        if (vctsize > 0)
          {
            double *vctptr = (double *) Malloc(vctsize * sizeof(double));
            size_t dummy = vctsize;
            GRIB_CHECK(grib_get_double_array(gh, "pv", vctptr, &dummy), 0);
            varDefVCT(vctsize, vctptr);
            Free(vctptr);
          }
        break;
      }
    case ZAXIS_REFERENCE:
      {
        unsigned char uuid[CDI_UUID_SIZE];
        long lpar;
        GRIB_CHECK(grib_get_long(gh, "NV", &lpar), 0);
        // if (lpar != 6) fprintf(stderr, "Warning ...\n");
        GRIB_CHECK(grib_get_long(gh, "nlev", &lpar), 0);
        int nhlev = (int) lpar;
        GRIB_CHECK(grib_get_long(gh, "numberOfVGridUsed", &lpar), 0);
        int nvgrid = (int) lpar;
        size_t len = (size_t) CDI_UUID_SIZE;
        memset(uuid, 0, CDI_UUID_SIZE);
        GRIB_CHECK(grib_get_bytes(gh, "uuidOfVGrid", uuid, &len), 0);
        varDefZAxisReference(nhlev, nvgrid, uuid);
        break;
      }
    }

  // if ( datatype > 32 ) datatype = CDI_DATATYPE_PACK32;
  if (datatype < 0) datatype = CDI_DATATYPE_PACK;

  // add the previously read record data to the (intermediate) list of records
  int tile_index = 0;
  int varID = 0, levelID = 0;
  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, level_sf, level_unit, datatype, &varID, &levelID,
               tsteptype, leveltype1, leveltype2, varname, scanKeys, tiles, &tile_index);

  record->varID = (short) varID;
  record->levelID = levelID;

  varDefCompType(varID, comptype);

  if (uvRelativeToGrid) varDefKeyInt(varID, CDI_KEY_UVRELATIVETOGRID, 1);

  if (varname[0]) gribapiGetNameKeys(gh, varID);
  gribapiGetKeys(gh, varID);

  if (lread_additional_keys)
    {
      long lval;
      double dval;
      for (int i = 0; i < cdiNAdditionalGRIBKeys; i++)
        {
          // note: if the key is not defined, we do not throw an error!
          if (grib_get_long(gh, cdiAdditionalGRIBKeys[i], &lval) == 0)
            varDefOptGribInt(varID, tile_index, lval, cdiAdditionalGRIBKeys[i]);
          if (grib_get_double(gh, cdiAdditionalGRIBKeys[i], &dval) == 0)
            varDefOptGribDbl(varID, tile_index, dval, cdiAdditionalGRIBKeys[i]);
        }
    }

  if (varInqInst(varID) == CDI_UNDEFID)
    {
      long center, subcenter;
      GRIB_CHECK(grib_get_long(gh, "centre", &center), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter), 0);
      int instID = institutInq((int) center, (int) subcenter, NULL, NULL);
      if (instID == CDI_UNDEFID) instID = institutDef((int) center, (int) subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if (varInqModel(varID) == CDI_UNDEFID)
    {
      long processID;
      if (grib_get_long(gh, "generatingProcessIdentifier", &processID) == 0)
        {
          /* FIXME: assert(processID >= INT_MIN && processID <= INT_MAX) */
          int modelID = modelInq(varInqInst(varID), (int) processID, NULL);
          if (modelID == CDI_UNDEFID) modelID = modelDef(varInqInst(varID), (int) processID, NULL);
          varDefModel(varID, modelID);
        }
    }

  if (varInqTable(varID) == CDI_UNDEFID)
    {
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);

      if (pdis == 255)
        {
          int tabnum = pcat;
          int tableID = tableInq(varInqModel(varID), tabnum, NULL);
          if (tableID == CDI_UNDEFID) tableID = tableDef(varInqModel(varID), tabnum, NULL);
          varDefTable(varID, tableID);
        }
    }

  if (CDI_Debug)
    Message("varID = %d  param = %d  zaxistype = %d  gridID = %d  levelID = %d", varID, param, zaxistype, gridID, levelID);
}

static compvar2_t
gribapiVarSet(int param, int level1, int level2, int leveltype, int tsteptype, size_t gridsize, char *name, VarScanKeys scanKeys,
              var_tile_t tiles_data)
{
  compvar2_t compVar;
  memset(&compVar, 0, sizeof(compvar2_t));
  size_t maxlen = sizeof(compVar.name);
  size_t len = strlen(name);
  if (len > maxlen) len = maxlen;

  compVar.param = param;
  compVar.level1 = level1;
  compVar.level2 = level2;
  compVar.ltype = leveltype;
  compVar.tsteptype = tsteptype;
  compVar.gridsize = gridsize;
  // memset(compVar.name, 0, maxlen);
  memcpy(compVar.name, name, len);
  compVar.scanKeys = scanKeys;
  compVar.tiles = tiles_data;

  return compVar;
}

static int
gribapiVarCompare(const compvar2_t *compVar, const record_t *record, int flag)
{
  compvar2_t compVar0;
  memset(&compVar0, 0, sizeof(compvar2_t));
  compVar0.param = record->param;
  compVar0.level1 = record->ilevel;
  compVar0.level2 = record->ilevel2;
  compVar0.ltype = record->ltype;
  compVar0.tsteptype = record->tsteptype;
  compVar0.gridsize = record->gridsize;
  memcpy(compVar0.name, record->varname, sizeof(compVar->name));

  if (flag == 0)
    {
      if (compVar0.tsteptype == TSTEP_INSTANT && compVar->tsteptype == TSTEP_INSTANT3) compVar0.tsteptype = TSTEP_INSTANT3;
      if (compVar0.tsteptype == TSTEP_INSTANT3 && compVar->tsteptype == TSTEP_INSTANT) compVar0.tsteptype = TSTEP_INSTANT;
    }

  compVar0.scanKeys = record->scanKeys;
  compVar0.tiles = record->tiles;

  // printf("var1: level1=%d level2=%d\n", compVar.level1, compVar.level2);
  // printf("var2: level1=%d level2=%d\n", compVar0.level1, compVar0.level2);

  return memcmp(&compVar0, compVar, sizeof(compvar2_t));
}

static grib_handle *
gribapiGetDiskRepresentation(size_t recsize, size_t *buffersize, void **gribbuffer, int *outDatatype, int *outCompressionType)
{
  int gribversion = (int) ((char *) *gribbuffer)[7];

  if (gribversion <= 1) *outCompressionType = grbDecompress(recsize, buffersize, gribbuffer);

  grib_handle *gh = grib_handle_new_from_message(NULL, *gribbuffer, recsize);

  bool lieee = false;

  if (gribversion > 1)
    {
      size_t len = 256;
      char typeOfPacking[256];
      if (grib_get_string(gh, "packingType", typeOfPacking, &len) == 0)
        {
          // fprintf(stderr, "packingType %zu %s\n", len, typeOfPacking);
          if (strncmp(typeOfPacking, "grid_jpeg", len) == 0)
            *outCompressionType = CDI_COMPRESS_JPEG;
          else if (strncmp(typeOfPacking, "grid_ccsds", len) == 0)
            *outCompressionType = CDI_COMPRESS_AEC;
          else if (strncmp(typeOfPacking, "grid_ieee", len) == 0)
            lieee = true;
        }
    }

  if (lieee)
    {
      long precision;
      int status = grib_get_long(gh, "precision", &precision);
      *outDatatype = (status == 0 && precision == 1) ? CDI_DATATYPE_FLT32 : CDI_DATATYPE_FLT64;
    }
  else
    {
      *outDatatype = CDI_DATATYPE_PACK;
      long bitsPerValue;
      if (grib_get_long(gh, "bitsPerValue", &bitsPerValue) == 0)
        {
          if (bitsPerValue > 0 && bitsPerValue <= 32) *outDatatype = (int) bitsPerValue;
        }
    }

  return gh;
}

typedef enum
{
  CHECKTIME_OK,
  CHECKTIME_SKIP,
  CHECKTIME_STOP,
  CHECKTIME_INCONSISTENT
} checkTimeResult;

static checkTimeResult
checkTime(stream_t *streamptr, const compvar2_t *compVar, CdiDateTime verificationTime, CdiDateTime expectedVTime)
{
  // First determine whether the current record exists already.
  int recID = 0;
  for (; recID < streamptr->nrecs; recID++)
    {
      if (gribapiVarCompare(compVar, &streamptr->tsteps[0].records[recID], 1) == 0) break;
    }
  bool recordExists = (recID < streamptr->nrecs);

  // Then we need to know whether the verification time is consistent.
  bool consistentTime = cdiDateTime_isEQ(verificationTime, expectedVTime);

  // Finally, we make a decision.
  if (CDI_Inventory_Mode == 1)
    {
      if (recordExists) return CHECKTIME_STOP;
      if (!consistentTime) return CHECKTIME_INCONSISTENT;
    }
  else
    {
      if (!consistentTime) return CHECKTIME_STOP;
      if (recordExists) return CHECKTIME_SKIP;
    }

  return CHECKTIME_OK;
}

#define gribWarning(text, nrecs, timestep, varname, param, level1, level2)                                                      \
  do                                                                                                                            \
    {                                                                                                                           \
      char paramstr[32];                                                                                                        \
      cdiParamToString(param, paramstr, sizeof(paramstr));                                                                      \
      Warning("Record %2d (name=%s id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, varname, paramstr, level1, level2, timestep, \
              text);                                                                                                            \
    }                                                                                                                           \
  while (0)

#ifdef HAVE_LIBFDB5
#include "cdi_fdb.h"
#endif

int
fdbScanTimesteps(stream_t *streamptr)
{
  (void) streamptr;
#ifdef HAVE_LIBFDB5
  void *gribbuffer = NULL;
  size_t buffersize = 0;
  grib_handle *gh = NULL;

  fdb_handle_t *fdbHandle = streamptr->protocolData;

  fdb_request_t *request = cdi_create_fdb_request(streamptr->filename);

  fdbKeyValueEntry *keyValueList = NULL;
  int numItems = cdi_fdb_fill_kvlist(fdbHandle, request, &keyValueList);
  fdb_delete_request(request);
  if (numItems == 0) Error("FDB request doesn't find any database entries!");
  if (CDI_Debug)
    {
      printf("Original FDB items:\n");
      print_keyvalueList(numItems, keyValueList);
    }

  // if (check_keyvalueList(numItems, keyValueList) != 0) return CDI_EUFSTRUCT;

  RecordInfoEntry *recordInfoList = (RecordInfoEntry *) malloc(numItems * sizeof(RecordInfoEntry));
  if (decode_keyvalue(numItems, keyValueList, recordInfoList) != 0) return CDI_EUFSTRUCT;

  cdi_fdb_sort_datetime(numItems, recordInfoList);

  if (CDI_Debug)
    {
      printf("Sorted FDB items:\n");
      print_keyvalueList_sorted(numItems, keyValueList, recordInfoList);
    }

  int numRecords = get_num_records(numItems, recordInfoList);
  if (numRecords == 0) return CDI_EUFSTRUCT;

  int numTimesteps = numItems / numRecords;
  if (CDI_Debug) Message("numRecords=%d  numTimesteps=%d", numRecords, numTimesteps);

  int *timestepRecordOffset = (int *) malloc(numTimesteps * sizeof(int));
  for (int i = 0; i < numTimesteps; i++) timestepRecordOffset[i] = i * numRecords;
  numTimesteps = remove_duplicate_timesteps(recordInfoList, numRecords, numTimesteps, timestepRecordOffset);
  if (CDI_Debug) Message("numRecords=%d  numTimesteps=%d", numRecords, numTimesteps);
  // Message("numRecords=%d  numTimesteps=%d", numRecords, numTimesteps);

  // CdiDateTime vDateTime0;
  // cdiDateTime_init(&vDateTime0);
  int fcast = 0;

  streamptr->curTsID = 0;

  int tsIDnew = tstepsNewEntry(streamptr);
  if (tsIDnew != 0) Error("Internal problem! tstepsNewEntry returns %d", tsIDnew);

  taxis_t *taxis = &streamptr->tsteps[tsIDnew].taxis;

  for (int recID = 0; recID < numRecords; recID++)
    {
      int fdbItem = recordInfoList[recID].fdbIndex;
      int vdate = recordInfoList[recID].date;
      int vtime = recordInfoList[recID].time * 100;

      long recsize = cdi_fdb_read_record(fdbHandle, &keyValueList[fdbItem], &buffersize, &gribbuffer);

      int datatype, comptype = 0;
      gh = gribapiGetDiskRepresentation(recsize, &buffersize, &gribbuffer, &datatype, &comptype);

      GRIB_CHECK(my_grib_set_double(gh, "missingValue", CDI_Default_Missval), 0);

      int leveltype1 = -1, leveltype2 = -1, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      int level1 = 0, level2 = 0;
      grib_get_level(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      gribapi_get_string(gh, "shortName", varname, sizeof(varname));

      int param = gribapiGetParam(gh);

      CdiDateTime sDateTime;
      // CdiDateTime vDateTime = gribapiGetValidityDateTime(gh, &sDateTime);
      CdiDateTime vDateTime = cdiDateTime_set(vdate, vtime);
      sDateTime = vDateTime;

      VarScanKeys scanKeys = gribapiGetScanKeys(gh);

      if (recID == 0)
        {
          // vDateTime0 = vDateTime;
          taxis->rDateTime = gribapiGetDataDateTime(gh);
          fcast = gribapiTimeIsFC(gh);
          if (fcast) taxis->unit = gribapi_get_timeunits(gh);
          taxis->fDateTime = taxis->rDateTime;
          taxis->sDateTime = sDateTime;
          taxis->vDateTime = vDateTime;
        }

      if (CDI_Debug)
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4d name=%s id=%s ltype=%d lev1=%d lev2=%d vDateTime=%s", recID + 1, varname, paramstr, leveltype1, level1,
                  level2, CdiDateTime_string(vDateTime));
        }

      var_tile_t *ptiles = memcmp(&tiles, &dummy_tiles, sizeof(var_tile_t)) ? &tiles : NULL;
      int recpos = 0;
      gribapiAddRecord(streamptr, param, gh, recsize, recpos, datatype, comptype, varname, leveltype1, leveltype2, lbounds, level1,
                       level2, level_sf, level_unit, &scanKeys, ptiles, true, fdbItem);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if (gh) grib_handle_delete(gh);

  cdi_generate_vars(streamptr);

  taxis->type = fcast ? TAXIS_RELATIVE : TAXIS_ABSOLUTE;
  int taxisID = taxisCreate(taxis->type);

  vlistDefTaxis(streamptr->vlistID, taxisID);

  streamScanResizeRecords1(streamptr);

  streamptr->record->buffer = gribbuffer;
  streamptr->record->buffersize = buffersize;

  if (numTimesteps == 1) streamptr->ntsteps = 1;
  streamScanTimeConstAdjust(streamptr, taxis);

  for (int tsID = 1; tsID < numTimesteps; tsID++)
    {
      int recordOffset = timestepRecordOffset[tsID];
      int vdate = recordInfoList[recordOffset].date;
      int vtime = recordInfoList[recordOffset].time * 100;
      // printf("timestep=%d recOffset=%d date=%d time=%d\n", tsID + 1, recordOffset, vdate, vtime);

      int tsIDnext = tstepsNewEntry(streamptr);
      if (tsIDnext != tsID) Error("Internal error. tsID = %d", tsID);

      streamptr->tsteps[tsID - 1].next = true;
      streamptr->tsteps[tsID].position = 0;

      taxis = &streamptr->tsteps[tsID].taxis;

      cdi_create_records(streamptr, tsID);
      record_t *records = streamptr->tsteps[tsID].records;

      int nrecs = (tsID == 1) ? streamScanInitRecords2(streamptr) : streamScanInitRecords(streamptr, tsID);
      if (nrecs != numRecords) Error("Internal error. nrecs = %d", nrecs);

      taxis->vDateTime = cdiDateTime_set(vdate, vtime);

      int rindex = 0;
      for (int recID = 0; recID < numRecords; recID++)
        {
          records[recID].used = true;
          streamptr->tsteps[tsID].recIDs[rindex] = recID;
          rindex++;

          records[recID].position = 0;
          records[recID].size = 0;
          records[recID].fdbItemIndex = recordInfoList[recordOffset + recID].fdbIndex;
        }

      if (tsID == 1) streamptr->tsteps[1].nrecs = numRecords;
    }

  streamptr->rtsteps = numTimesteps;
  streamptr->ntsteps = numTimesteps;

  streamptr->fdbNumItems = numItems;
  streamptr->fdbKeyValueList = keyValueList;

  if (recordInfoList) free(recordInfoList);
  if (timestepRecordOffset) free(timestepRecordOffset);
#endif

  return 0;
}

int
gribapiScanTimestep1(stream_t *streamptr)
{
  CdiDateTime vDateTime0;
  cdiDateTime_init(&vDateTime0);
  off_t recpos = 0;
  void *gribbuffer = NULL;
  size_t buffersize = 0;
  int nrecsScanned = 0;  // Only used for debug output.
  bool warn_time = true;
  int fcast = 0;
  grib_handle *gh = NULL;

  streamptr->curTsID = 0;

  int tsID = tstepsNewEntry(streamptr);
  if (tsID != 0) Error("Internal problem! tstepsNewEntry returns %d", tsID);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  int fileID = streamptr->fileID;

  unsigned nrecs = 0;
  while (true)
    {
      size_t recsize = gribGetSize(fileID);
      recpos = fileGetPos(fileID);

      if (recsize == 0)
        {
          streamptr->ntsteps = 1;
          break;
        }

      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      // Search for next 'GRIB', read the following record, and position file offset after it.
      if (gribRead(fileID, gribbuffer, &readsize)) break;

      nrecsScanned++;

      int datatype, comptype = 0;
      gh = gribapiGetDiskRepresentation(recsize, &buffersize, &gribbuffer, &datatype, &comptype);

      GRIB_CHECK(my_grib_set_double(gh, "missingValue", CDI_Default_Missval), 0);

      int leveltype1 = -1, leveltype2 = -1, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      int level1 = 0, level2 = 0;
      grib_get_level(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      gribapi_get_string(gh, "shortName", varname, sizeof(varname));

      int param = gribapiGetParam(gh);

      if (!varname[0]) param_to_name(param, varname);

      CdiDateTime sDateTime;
      CdiDateTime vDateTime = gribapiGetValidityDateTime(gh, &sDateTime);

      VarScanKeys scanKeys = gribapiGetScanKeys(gh);

      if (nrecs == 0)
        {
          vDateTime0 = vDateTime;
          taxis->rDateTime = gribapiGetDataDateTime(gh);
          fcast = gribapiTimeIsFC(gh);
          if (fcast) taxis->unit = gribapi_get_timeunits(gh);
          taxis->fDateTime = taxis->rDateTime;
          taxis->sDateTime = sDateTime;
          taxis->vDateTime = vDateTime;
        }
      else
        {
          if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;

          int tsteptype = gribapiGetTsteptype(gh);
          size_t gridsize = gribapiGetGridsize(gh);
          compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, gridsize, varname, scanKeys, tiles);
          checkTimeResult result = checkTime(streamptr, &compVar, vDateTime, vDateTime0);
          if (result == CHECKTIME_STOP)
            {
              nrecsScanned--;
              break;
            }
          else if (result == CHECKTIME_SKIP)
            {
              gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, varname, param, level1, level2);
              continue;
            }
          else if (result == CHECKTIME_INCONSISTENT && warn_time)
            {
              gribWarning("Inconsistent verification time!", nrecsScanned, tsID + 1, varname, param, level1, level2);
              warn_time = false;
            }
          assert(result == CHECKTIME_OK || result == CHECKTIME_INCONSISTENT);
        }

      nrecs++;

      if (CDI_Debug)
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4u %8d name=%s id=%s ltype=%d lev1=%d lev2=%d vDateTime=%s", nrecs, (int) recpos, varname, paramstr, leveltype1,
                  level1, level2, CdiDateTime_string(vDateTime));
        }

      var_tile_t *ptiles = memcmp(&tiles, &dummy_tiles, sizeof(var_tile_t)) ? &tiles : NULL;
      gribapiAddRecord(streamptr, param, gh, recsize, recpos, datatype, comptype, varname, leveltype1, leveltype2, lbounds, level1,
                       level2, level_sf, level_unit, &scanKeys, ptiles, true, -1);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if (gh) grib_handle_delete(gh);

  streamptr->rtsteps = 1;

  if (nrecs == 0) return CDI_EUFSTRUCT;

  if (streamptr->query)
    {
      int numEntries = cdiQueryNumEntries(streamptr->query);
      int numEntriesFound = cdiQueryNumEntriesFound(streamptr->query);
      cdiQueryPrintEntriesNotFound(streamptr->query);
      if (numEntriesFound == 0 || (CDI_Query_Abort && numEntries != numEntriesFound)) return CDI_EQENF;
    }

  cdi_generate_vars(streamptr);

  taxis->type = fcast ? TAXIS_RELATIVE : TAXIS_ABSOLUTE;
  int taxisID = taxisCreate(taxis->type);

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  streamScanResizeRecords1(streamptr);

  streamptr->record->buffer = gribbuffer;
  streamptr->record->buffersize = buffersize;

  streamScanTsFixNtsteps(streamptr, recpos);
  streamScanTimeConstAdjust(streamptr, taxis);

  return 0;
}

int
gribapiScanTimestep2(stream_t *streamptr)
{
  CdiDateTime vDateTime0;
  cdiDateTime_init(&vDateTime0);
  int rstatus = 0;
  off_t recpos = 0;
  // int gridID;
  int recID;
  grib_handle *gh = NULL;

  streamptr->curTsID = 1;

  int fileID = streamptr->fileID;
  int vlistID = streamptr->vlistID;

  void *gribbuffer = streamptr->record->buffer;
  size_t buffersize = streamptr->record->buffersize;

  int tsID = streamptr->rtsteps;
  if (tsID != 1) Error("Internal problem! unexpected timestep %d", tsID + 1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);
  record_t *records = streamptr->tsteps[tsID].records;

  int nrecords = streamScanInitRecords2(streamptr);

  int nrecsScanned = nrecords;  // Only used for debug output
  for (int rindex = 0; rindex <= nrecords; ++rindex)
    {
      size_t recsize = gribGetSize(fileID);
      recpos = fileGetPos(fileID);
      if (recsize == 0)
        {
          streamptr->ntsteps = 2;
          break;
        }

      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      if (gribRead(fileID, gribbuffer, &readsize)) break;

      grbDecompress(recsize, &buffersize, &gribbuffer);

      nrecsScanned++;
      gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", CDI_Default_Missval), 0);

      int level1 = 0, level2 = 0, leveltype1, leveltype2, lbounds, level_sf, level_unit;
      var_tile_t tiles = dummy_tiles;
      grib_get_level(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

      char varname[256];
      gribapi_get_string(gh, "shortName", varname, sizeof(varname));

      int param = gribapiGetParam(gh);

      if (!varname[0]) param_to_name(param, varname);

      CdiDateTime sDateTime;
      CdiDateTime vDateTime = gribapiGetValidityDateTime(gh, &sDateTime);

      if (rindex == 0)
        {
          vDateTime0 = vDateTime;
          int taxisID = vlistInqTaxis(vlistID);
          if (taxisInqType(taxisID) == TAXIS_RELATIVE)
            {
              taxis->type = TAXIS_RELATIVE;
              taxis->unit = gribapi_get_timeunits(gh);
              taxis->rDateTime = gribapiGetDataDateTime(gh);
            }
          else
            {
              taxis->type = TAXIS_ABSOLUTE;
            }
          taxis->fDateTime = taxis->rDateTime;
          taxis->vDateTime = vDateTime;
          taxis->sDateTime = sDateTime;
        }
      else
        {
          if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;
        }

      VarScanKeys scanKeys = gribapiGetScanKeys(gh);

      int tsteptype = gribapiGetTsteptype(gh);
      size_t gridsize = gribapiGetGridsize(gh);
      compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, gridsize, varname, scanKeys, tiles);

      for (recID = 0; recID < nrecords; recID++)
        if (gribapiVarCompare(&compVar, &records[recID], 0) == 0) break;

      if (recID == nrecords)
        {
          if (CDI_Inventory_Mode == 1)
            {
              gribWarning("Parameter not defined at timestep 1!", nrecsScanned, tsID + 1, varname, param, level1, level2);
              return CDI_EUFSTRUCT;
            }
          else
            {
              gribWarning("Parameter not defined at timestep 1, skipped!", nrecsScanned, tsID + 1, varname, param, level1, level2);
              continue;
            }
        }

      if (records[recID].used)
        {
          if (CDI_Inventory_Mode == 1)
            break;
          else
            {
              if (cdiDateTime_isNE(vDateTime, vDateTime0)) break;

              gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, varname, param, level1, level2);
              continue;
            }
        }

      records[recID].used = true;
      streamptr->tsteps[tsID].recIDs[rindex] = recID;

      if (CDI_Debug)
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Message("%4d %8lld name=%s id=%s ltype=%d lev1=%d lev2=%d vDateTime=%s", nrecsScanned, (long long) recpos, varname,
                  paramstr, leveltype1, level1, level2, CdiDateTime_string(vDateTime));
        }

      if (gribapiVarCompare(&compVar, &records[recID], 0) != 0)
        {
          Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                  records[recID].ilevel, level1);
          return CDI_EUFSTRUCT;
        }

      records[recID].position = recpos;
      records[recID].size = recsize;

      int varID = records[recID].varID;

      if (tsteptype != vlistInqVarTsteptype(vlistID, varID)) vlistDefVarTsteptype(vlistID, varID, tsteptype);

      grib_handle_delete(gh);
      gh = NULL;
    }

  if (gh) grib_handle_delete(gh);

  int nrecs = 0;
  for (recID = 0; recID < nrecords; recID++)
    {
      if (records[recID].used)
        nrecs++;
      else
        vlistDefVarTimetype(vlistID, records[recID].varID, TIME_CONSTANT);
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  streamScanTsFixNtsteps(streamptr, recpos);

  streamptr->record->buffer = gribbuffer;
  streamptr->record->buffersize = buffersize;

  return rstatus;
}

int
gribapiScanTimestep(stream_t *streamptr)
{
  int vrecID, recID = -1;
  int nrecs = 0;
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  if (streamptr->tsteps[tsID].recordSize == 0)
    {
      void *gribbuffer = streamptr->record->buffer;
      size_t buffersize = streamptr->record->buffersize;

      cdi_create_records(streamptr, tsID);
      record_t *records = streamptr->tsteps[tsID].records;

      nrecs = streamScanInitRecords(streamptr, tsID);

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      int nrecsScanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs * (tsID - 1);  // Only used for debug output.
      off_t recpos = 0;
      CdiDateTime vDateTime0;
      cdiDateTime_init(&vDateTime0);
      grib_handle *gh = NULL;
      char varname[256];
      for (int rindex = 0; rindex <= nrecs; ++rindex)
        {
          varname[0] = 0;
          size_t recsize = gribGetSize(fileID);
          recpos = fileGetPos(fileID);
          if (recsize == 0)
            {
              streamptr->ntsteps = streamptr->rtsteps + 1;
              break;
            }

          ensureBufferSize(recsize, &buffersize, &gribbuffer);

          size_t readsize = recsize;
          if (gribRead(fileID, gribbuffer, &readsize))
            {
              Warning("Inconsistent timestep %d (GRIB record %d/%d)!", tsID + 1, rindex + 1, streamptr->tsteps[tsID].recordSize);
              break;
            }

          grbDecompress(recsize, &buffersize, &gribbuffer);

          nrecsScanned++;
          gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
          GRIB_CHECK(my_grib_set_double(gh, "missingValue", CDI_Default_Missval), 0);

          int level1 = 0, level2 = 0, leveltype1, leveltype2 = -1, lbounds, level_sf, level_unit;
          var_tile_t tiles = dummy_tiles;
          grib_get_level(gh, &leveltype1, &leveltype2, &lbounds, &level1, &level2, &level_sf, &level_unit, &tiles);

          CdiDateTime sDateTime;
          CdiDateTime vDateTime = gribapiGetValidityDateTime(gh, &sDateTime);

          if (rindex == nrecs) break;

          gribapi_get_string(gh, "shortName", varname, sizeof(varname));

          int param = gribapiGetParam(gh);

          if (!varname[0]) param_to_name(param, varname);

          if (rindex == 0)
            {
              vDateTime0 = vDateTime;
              int taxisID = vlistInqTaxis(vlistID);
              if (taxisInqType(taxisID) == TAXIS_RELATIVE)
                {
                  taxis->type = TAXIS_RELATIVE;
                  taxis->unit = gribapi_get_timeunits(gh);
                  taxis->rDateTime = gribapiGetDataDateTime(gh);
                }
              else
                {
                  taxis->type = TAXIS_ABSOLUTE;
                }
              taxis->fDateTime = taxis->rDateTime;
              taxis->vDateTime = vDateTime;
              taxis->sDateTime = sDateTime;
            }
          else
            {
              if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;
            }

          VarScanKeys scanKeys = gribapiGetScanKeys(gh);

          int tsteptype = gribapiGetTsteptype(gh);
          size_t gridsize = gribapiGetGridsize(gh);
          compvar2_t compVar = gribapiVarSet(param, level1, level2, leveltype1, tsteptype, gridsize, varname, scanKeys, tiles);

          for (vrecID = 0; vrecID < nrecs; vrecID++)
            {
              recID = streamptr->tsteps[1].recIDs[vrecID];
              if (gribapiVarCompare(&compVar, &records[recID], 0) == 0) break;
            }

          if (vrecID == nrecs)
            {
              if (CDI_Inventory_Mode == 1)
                {
                  gribWarning("Parameter not defined at timestep 1!", nrecsScanned, tsID + 1, varname, param, level1, level2);
                  return CDI_EUFSTRUCT;
                }
              else
                {
                  gribWarning("Parameter not defined at timestep 1, skipped!", nrecsScanned, tsID + 1, varname, param, level1,
                              level2);
                  continue;
                }
            }

          if (CDI_Inventory_Mode != 1)
            {
              if (records[recID].used)
                {
                  if (cdiDateTime_isNE(vDateTime, vDateTime0)) break;

                  if (CDI_Debug)
                    gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, varname, param, level1, level2);

                  continue;
                }
            }

          records[recID].used = true;
          streamptr->tsteps[tsID].recIDs[rindex] = recID;

          if (CDI_Debug)
            Message("%4d %8lld %4d %8d %8s", rindex + 1, (long long) recpos, param, level1, CdiDateTime_string(vDateTime));

          if (gribapiVarCompare(&compVar, &records[recID], 0) != 0)
            {
              Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                      records[recID].ilevel, level1);
              Error("Invalid, unsupported or inconsistent record structure");
            }

          records[recID].position = recpos;
          records[recID].size = recsize;

          if (CDI_Debug) Message("%4d %8lld %4d %8d %s", rindex, (long long) recpos, param, level1, CdiDateTime_string(vDateTime));

          grib_handle_delete(gh);
          gh = NULL;
        }

      if (gh) grib_handle_delete(gh);

      for (vrecID = 0; vrecID < nrecs; vrecID++)
        {
          recID = streamptr->tsteps[tsID].recIDs[vrecID];
          if (!records[recID].used) break;
        }

      if (vrecID < nrecs)
        {
          gribWarning("Parameter not found!", nrecsScanned, tsID + 1, varname, records[recID].param, records[recID].ilevel,
                      records[recID].ilevel2);
          return CDI_EUFSTRUCT;
        }

      streamptr->rtsteps++;

      if (streamptr->ntsteps != streamptr->rtsteps)
        {
          tsID = tstepsNewEntry(streamptr);
          if (tsID != streamptr->rtsteps) Error("Internal error. tsID = %d", tsID);

          streamptr->tsteps[tsID - 1].next = true;
          streamptr->tsteps[tsID].position = recpos;
        }

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);
      streamptr->tsteps[tsID].position = recpos;

      streamptr->record->buffer = gribbuffer;
      streamptr->record->buffersize = buffersize;
    }

  if (nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs)
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return streamptr->ntsteps;
}

#ifdef gribWarning
#undef gribWarning
#endif

static void
unpack_alternative_rows(grib_handle *gh, int memType, void *data)
{
  long xsize = 0, ysize = 0;
  grib_get_long(gh, "Nx", &xsize);
  grib_get_long(gh, "Ny", &ysize);

  if (memType == MEMTYPE_FLOAT)
    {
      float *pdata = (float *) data;
      for (int j = 1; j < ysize; j += 2)
        for (int i = 0; i < xsize / 2; i++)
          {
            float tmp = pdata[j * xsize + i];
            pdata[j * xsize + i] = pdata[j * xsize + xsize - i - 1];
            pdata[j * xsize + xsize - i - 1] = tmp;
          }
    }
  else
    {
      double *pdata = (double *) data;
      for (int j = 1; j < ysize; j += 2)
        for (int i = 0; i < xsize / 2; i++)
          {
            double tmp = pdata[j * xsize + i];
            pdata[j * xsize + i] = pdata[j * xsize + xsize - i - 1];
            pdata[j * xsize + xsize - i - 1] = tmp;
          }
    }
}

int
gribapiDecode(int memType, void *gribbuffer, size_t gribsize, void *data, size_t gridsize, int unreduced, size_t *numMissVals,
              double missval)
{
  int status = 0;

  if (unreduced)
    {
      static bool lwarn = true;
      if (lwarn)
        {
          lwarn = false;
          Warning("Conversion of gaussian reduced grids unsupported!");
        }
    }

  size_t recsize = (size_t) gribsize;
  grib_handle *gh = grib_handle_new_from_message(NULL, gribbuffer, recsize);
  GRIB_CHECK(my_grib_set_double(gh, "missingValue", missval), 0);

  // get the size of the values array
  size_t datasize;
  GRIB_CHECK(grib_get_size(gh, "values", &datasize), 0);
  // long numberOfPoints;
  // GRIB_CHECK(grib_get_long(gh, "numberOfPoints", &numberOfPoints), 0);
  // printf("values_size = %d  numberOfPoints = %ld\n", datasize, numberOfPoints);

  if (datasize != gridsize) Error("numberOfPoint (%zu) and gridSize (%zu) differ!", datasize, gridsize);
  size_t dummy = datasize;

  if (memType == MEMTYPE_FLOAT)
    {
#ifdef HAVE_GRIBAPI_FLOAT_INTERFACE
      GRIB_CHECK(grib_get_float_array(gh, "values", (float *) data, &dummy), 0);
#else
      Error("grib_get_float_array() not found!");
#endif
    }
  else
    {
      GRIB_CHECK(grib_get_double_array(gh, "values", (double *) data, &dummy), 0);
    }

  if (gribEditionNumber(gh) > 1)
    {
      long alternativeRowScanning = false;
      grib_get_long(gh, "alternativeRowScanning", &alternativeRowScanning);
      if (alternativeRowScanning) unpack_alternative_rows(gh, memType, data);
    }

  long lpar;
  GRIB_CHECK(grib_get_long(gh, "gridDefinitionTemplateNumber", &lpar), 0);
  int gridtype = (int) lpar;

  *numMissVals = 0;
  if (gridtype < 50 || gridtype > 53)
    {
      GRIB_CHECK(grib_get_long(gh, "numberOfMissing", &lpar), 0);
      *numMissVals = (int) lpar;
      // printf("gridtype %d, numMissVals %d\n", gridtype, numMissVals);
    }

  grib_handle_delete(gh);

  return status;
}

static void
gribapiDefInstitut(grib_handle *gh, int vlistID, int varID)
{
  int instID = vlistInqInstitut(vlistID);

  if (instID == CDI_UNDEFID) instID = vlistInqVarInstitut(vlistID, varID);

  if (instID != CDI_UNDEFID)
    {
      long center = institutInqCenter(instID);
      long subcenter = institutInqSubcenter(instID);

      long center0, subcenter0;
      GRIB_CHECK(grib_get_long(gh, "centre", &center0), 0);
      GRIB_CHECK(grib_get_long(gh, "subCentre", &subcenter0), 0);

      if (center != center0) GRIB_CHECK(my_grib_set_long(gh, "centre", center), 0);
      if (subcenter != subcenter0) GRIB_CHECK(my_grib_set_long(gh, "subCentre", subcenter), 0);
    }

  int status;
  int centre, subCentre;
  status = cdiInqKeyInt(vlistID, CDI_GLOBAL, CDI_KEY_CENTRE, &centre);
  if (status == 0) grib_set_long(gh, "centre", centre);
  status = cdiInqKeyInt(vlistID, CDI_GLOBAL, CDI_KEY_SUBCENTRE, &subCentre);
  if (status == 0) grib_set_long(gh, "subCentre", subCentre);

  status = cdiInqKeyInt(vlistID, varID, CDI_KEY_CENTRE, &centre);
  if (status == 0) grib_set_long(gh, "centre", centre);
  status = cdiInqKeyInt(vlistID, varID, CDI_KEY_SUBCENTRE, &subCentre);
  if (status == 0) grib_set_long(gh, "subCentre", subCentre);
}

static void
gribapiDefModel(grib_handle *gh, int vlistID, int varID)
{
  int modelID = vlistInqModel(vlistID);
  if (modelID == CDI_UNDEFID) modelID = vlistInqVarModel(vlistID, varID);

  if (modelID != CDI_UNDEFID) GRIB_CHECK(my_grib_set_long(gh, "generatingProcessIdentifier", modelInqGribID(modelID)), 0);
}

static void
gribapiDefParam(int editionNumber, grib_handle *gh, int param, const char *name, const char *stdname)
{
  bool ldefined = false;

  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  if (pnum < 0)
    {
      size_t len = strlen(stdname);
      if (len)
        {
          int status = my_grib_set_string(gh, "cfName", stdname, &len);
          if (status == 0)
            ldefined = true;
          else
            Warning("grib_api: No match for cfName=%s", stdname);
        }

      if (ldefined == false)
        {
          len = strlen(name);
          int status = my_grib_set_string(gh, "shortName", name, &len);
          if (status == 0)
            ldefined = true;
          else
            Warning("grib_api: No match for shortName=%s", name);
        }
    }

  if (ldefined == false)
    {
      if (pnum < 0) pnum = -pnum;

      if (pnum > 255)
        {
          static bool lwarn_pnum = true;
          if (lwarn_pnum)
            {
              Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum % 256);
              lwarn_pnum = false;
            }
          pnum = pnum % 256;
        }

      if (editionNumber <= 1)
        {
          static bool lwarn_pdis = true;
          if (pdis != 255 && lwarn_pdis)
            {
              char paramstr[32];
              cdiParamToString(param, paramstr, sizeof(paramstr));
              Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
              lwarn_pdis = false;
            }

          GRIB_CHECK(my_grib_set_long(gh, "table2Version", pcat), 0);
          GRIB_CHECK(my_grib_set_long(gh, "indicatorOfParameter", pnum), 0);
        }
      else
        {
          GRIB_CHECK(my_grib_set_long(gh, "discipline", pdis), 0);
          GRIB_CHECK(my_grib_set_long(gh, "parameterCategory", pcat), 0);
          GRIB_CHECK(my_grib_set_long(gh, "parameterNumber", pnum), 0);
        }
    }

  // printf("param: %d.%d.%d %s\n", pnum, pcat, pdis, name);
}

static int
getTimeunitFactor(int timeunit)
{
  switch (timeunit)
    {
    case TUNIT_SECOND: return 1;
    case TUNIT_MINUTE: return 60;
    case TUNIT_HOUR: return 3600;
    case TUNIT_3HOURS: return 10800;
    case TUNIT_6HOURS: return 21600;
    case TUNIT_12HOURS: return 43200;
    case TUNIT_DAY: return 86400;
    }

  return 3600;
}

static int
grib2ProDefTempHasStatisticalDef(int proDefTempNum)
{
  switch (proDefTempNum)
    {
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
    case 40034: return 1;
    }

  return 0;
}

static int
getUnitsOfTime(int timeunit)
{
  switch (timeunit)
    {
    case TUNIT_SECOND: return 13;
    case TUNIT_MINUTE: return 0;
    case TUNIT_HOUR: return 1;
    case TUNIT_3HOURS: return 10;
    case TUNIT_6HOURS: return 11;
    case TUNIT_12HOURS: return 12;
    case TUNIT_DAY: return 2;
    }

  return 1;
}

static void
gribapiDefStepUnits(int editionNumber, grib_handle *gh, int timeunit, int proDefTempNum, int gcinit)
{
  if (!gcinit)
    {
      long unitsOfTime = getUnitsOfTime(timeunit);

      grib_set_long(gh, "stepUnits", unitsOfTime);
      if (editionNumber == 1)
        {
          grib_set_long(gh, "unitOfTimeRange", unitsOfTime);
        }
      else if (grib2ProDefTempHasStatisticalDef(proDefTempNum))
        {
          grib_set_long(gh, "indicatorOfUnitForTimeRange", unitsOfTime);
          grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime);
        }
      else
        {
          // NOTE KNMI:  HIRLAM model files LAMH_D11 are in grib1 and do NOT have key indicatorOfUnitForTimeRange
          // Watch out for compatibility issues.
          grib_set_long(gh, "indicatorOfUnitOfTimeRange", unitsOfTime);
        }
    }
}

static int
gribapiDefSteptype(int editionNumber, grib_handle *gh, int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype,
                   int gcinit)
{
  const char *stepType = "instant";
  long proDefTempNum = 0;

  if (tsteptype >= TSTEP_INSTANT && tsteptype <= TSTEP_SUM)
    {
      stepType = cdiGribAPI_ts_str_map[tsteptype].sname;
      proDefTempNum = cdiGribAPI_ts_str_map[tsteptype].productionTemplate;
    }

  if (productDefinitionTemplate != -1)
    proDefTempNum = productDefinitionTemplate;
  else if (typeOfGeneratingProcess == 4)
    proDefTempNum = (proDefTempNum == 8) ? 11 : 1;

  if (!gcinit)
    {
      if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "productDefinitionTemplateNumber", proDefTempNum), 0);
      size_t len = strlen(stepType);
      int status = my_grib_set_string(gh, "stepType", stepType, &len);
      if (status != 0) GRIB_CHECK(my_grib_set_long(gh, "productDefinitionTemplateNumber", 0), 0);
    }

  return (int) proDefTempNum;
}

static void
gribapiDefDateTimeAbs(int editionNumber, grib_handle *gh, CdiDateTime dateTime, int productDefinitionTemplate,
                      int typeOfGeneratingProcess, int tsteptype, int gcinit)
{
  (void) gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);

  if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 0), 0);
  if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

  if (cdiDateTime_isNull(dateTime)) dateTime.date = cdiDate_set(10101);
  gribapiSetDataDateTime(gh, dateTime);
}

static int
gribapiDefDateTimeRel(int editionNumber, grib_handle *gh, CdiDateTime fDateTime, CdiDateTime vDateTime, CdiDateTime sDateTime,
                      int productDefinitionTemplate, int typeOfGeneratingProcess, int tsteptype, int timeunit, int calendar,
                      int gcinit)
{
  int status = -1;

  JulianDate julianDate1 = julianDate_encode(calendar, fDateTime);

  if (cdiDateTime_isNull(vDateTime)) vDateTime = fDateTime;

  JulianDate julianDate2 = julianDate_encode(calendar, vDateTime);
  JulianDate julianDate = julianDate_sub(julianDate2, julianDate1);

  int factor = getTimeunitFactor(timeunit);

  if (!(int) (fmod(julianDate_to_seconds(julianDate), factor)))
    {
      int proDefTempNum
          = gribapiDefSteptype(editionNumber, gh, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);
      gribapiDefStepUnits(editionNumber, gh, timeunit, proDefTempNum, gcinit);

      long startStep = 0;
      double endStepF = julianDate_to_seconds(julianDate) / factor;
      long maxStep = (editionNumber > 1) ? INT_MAX : 65000;
      if (endStepF > maxStep) return status;
      long endStep = lround(endStepF);

      bool hasStartDate = (tsteptype == TSTEP_RANGE || tsteptype == TSTEP_AVG || tsteptype == TSTEP_ACCUM || tsteptype == TSTEP_DIFF
                           || tsteptype == TSTEP_MIN || tsteptype == TSTEP_MAX || tsteptype == TSTEP_RMS || tsteptype == TSTEP_SD
                           || tsteptype == TSTEP_COV || tsteptype == TSTEP_RATIO || tsteptype == TSTEP_SUM);
      if (!cdiDateTime_isNull(sDateTime) && hasStartDate)
        {
          julianDate2 = julianDate_encode(calendar, sDateTime);
          startStep = lround(julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / factor);
        }

      if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "significanceOfReferenceTime", 1), 0);
      if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "stepRange", 0), 0);

      if (cdiDateTime_isNull(fDateTime)) fDateTime.date = cdiDate_set(10101);
      gribapiSetDataDateTime(gh, fDateTime);

      // printf(">>>>> tsteptype %d  startStep %ld  endStep %ld\n", tsteptype, startStep, endStep);

      // Product Definition Template Number: defined in GRIB_API file 4.0.table point in time products:
      if ((proDefTempNum >= 0 && proDefTempNum <= 7) || proDefTempNum == 55 || proDefTempNum == 40055)  // Tile
        startStep = endStep;

      if (endStep < startStep) return status;

      if (editionNumber == 1 && (startStep > 255 || endStep > 255))
        {
          startStep = 0;
          endStep = 0;
        }
      else
        {
          status = 0;
        }

      if (editionNumber > 1) GRIB_CHECK(my_grib_set_long(gh, "forecastTime", startStep), 0);
      // if ( editionNumber == 1 && startStep > 0) GRIB_CHECK(my_grib_set_long(gh, "startStep", startStep), 0);
      if (editionNumber == 1) GRIB_CHECK(my_grib_set_long(gh, "startStep", startStep), 0);
      GRIB_CHECK(my_grib_set_long(gh, "endStep", endStep), 0);
    }

  return status;
}

static void
gribapiDefTime(int editionNumber, int productDefinitionTemplate, int typeOfGeneratingProcess, grib_handle *gh,
               CdiDateTime vDateTime, int tsteptype, int numavg, int taxisID, int gcinit)
{
  UNUSED(numavg);

  int taxistype = (taxisID == -1) ? TAXIS_ABSOLUTE : taxisInqType(taxisID);

  if (typeOfGeneratingProcess == 196)
    {
      vDateTime = cdiDateTime_set(10101, 0);
      taxistype = TAXIS_ABSOLUTE;
    }

  if (taxistype == TAXIS_RELATIVE)
    {
      int timeunit = taxisInqTunit(taxisID);
      int calendar = taxisInqCalendar(taxisID);

      CdiDateTime fDateTime = taxisInqFdatetime(taxisID);
      if (cdiDateTime_isNull(fDateTime)) fDateTime = taxisInqRdatetime(taxisID);
      if (cdiDateTime_isLT(vDateTime, fDateTime)) fDateTime = vDateTime;

      CdiDateTime sDateTime = taxisInqSdatetime(taxisID);

      int status = gribapiDefDateTimeRel(editionNumber, gh, fDateTime, vDateTime, sDateTime, productDefinitionTemplate,
                                         typeOfGeneratingProcess, tsteptype, timeunit, calendar, gcinit);
      if (status != 0) taxistype = TAXIS_ABSOLUTE;
    }

  if (taxistype == TAXIS_ABSOLUTE)
    gribapiDefDateTimeAbs(editionNumber, gh, vDateTime, productDefinitionTemplate, typeOfGeneratingProcess, tsteptype, gcinit);
}

static void
gribapiDefGridRegular(grib_handle *gh, int gridID, int gridtype, bool gridIsRotated, bool gridIsCurvilinear, int uvRelativeToGrid)
{
  const char *mesg;
  // clang-format off
  if      (gridtype == GRID_GAUSSIAN)         mesg = "regular_gg";
  else if (gridtype == GRID_GAUSSIAN_REDUCED) mesg = "reduced_gg";
  else if (gridIsRotated)                     mesg = "rotated_ll";
  else                                        mesg = "regular_ll";
  // clang-format on
  size_t len = strlen(mesg);
  GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

  double xfirst = 0.0, xlast = 0.0, xinc = 0.0;
  double yfirst = 0.0, ylast = 0.0, yinc = 0.0;

  size_t nlon = gridInqXsize(gridID);
  size_t nlat = gridInqYsize(gridID);

  if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      xfirst = (nlon == 2) ? gridInqXval(gridID, 0) : 0.0;
      xlast = (nlon == 2) ? gridInqXval(gridID, 1) : 360.0 - 360.0 * 0.5 / (double) nlat;

      nlon = 0;

      int *reducedPoints = (int *) Malloc(nlat * sizeof(int));
      long *pl = (long *) Malloc(nlat * sizeof(long));
      gridInqReducedPoints(gridID, reducedPoints);
      for (size_t i = 0; i < nlat; ++i) pl[i] = reducedPoints[i];

      GRIB_CHECK(grib_set_long_array(gh, "pl", pl, nlat), 0);

      Free(pl);
      Free(reducedPoints);
    }
  else
    {
      if (nlon == 0)
        nlon = 1;
      else
        {
          xfirst = gridInqXval(gridID, 0);
          xlast = gridInqXval(gridID, (gridIsCurvilinear ? nlon * nlat : nlon) - 1);
          xinc = fabs(gridInqXinc(gridID));
        }
    }

  if (nlat == 0)
    nlat = 1;
  else
    {
      yfirst = gridInqYval(gridID, 0);
      ylast = gridInqYval(gridID, (gridIsCurvilinear ? nlon * nlat : nlat) - 1);
      yinc = fabs(gridInqYinc(gridID));
    }

  double xfirsto = xfirst;
  double xlasto = xlast;
  while (xfirsto > 360.0) xfirsto -= 360.0;
  while (xlasto > 360.0) xlasto -= 360.0;

  if (gridtype != GRID_GAUSSIAN_REDUCED) GRIB_CHECK(my_grib_set_long(gh, "Ni", nlon), 0);
  GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", xfirsto), 0);
  GRIB_CHECK(my_grib_set_double(gh, "longitudeOfLastGridPointInDegrees", xlasto), 0);
  if (gridtype != GRID_GAUSSIAN_REDUCED) GRIB_CHECK(my_grib_set_double(gh, "iDirectionIncrementInDegrees", xinc), 0);

  GRIB_CHECK(my_grib_set_long(gh, "Nj", (long) nlat), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", yfirst), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfLastGridPointInDegrees", ylast), 0);

  if (uvRelativeToGrid >= 0) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);

  GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", (xfirst > xlast)), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", (yfirst < ylast)), 0);

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      int np = gridInqNP(gridID);
      if (np == 0) np = nlat / 2;
      GRIB_CHECK(my_grib_set_long(gh, "numberOfParallelsBetweenAPoleAndTheEquator", np), 0);
    }
  else
    {
      GRIB_CHECK(my_grib_set_double(gh, "jDirectionIncrementInDegrees", yinc), 0);
    }

  if (gridIsRotated)
    {
      double xpole = 0.0, ypole = 0.0, angle = 0.0;
      gridInqParamRLL(gridID, &xpole, &ypole, &angle);

      xpole += 180.0;
      if (fabs(ypole) > 0.0) ypole = -ypole;  // change from north to south pole
      if (fabs(angle) > 0.0) angle = -angle;
      GRIB_CHECK(my_grib_set_double(gh, "latitudeOfSouthernPoleInDegrees", ypole), 0);
      GRIB_CHECK(my_grib_set_double(gh, "longitudeOfSouthernPoleInDegrees", xpole), 0);
      GRIB_CHECK(my_grib_set_double(gh, "angleOfRotation", angle), 0);
    }
}

static int
encode_shapeOfTheEarth(struct CDI_GridProjParams *gpp)
{
  int shapeOfTheEarth = 1;
  int a = (int) lround(gpp->a);
  int b = (int) lround(gpp->b);
  int rf = (int) lround(gpp->rf);

  // clang-format off
  if      (a == 6367470) shapeOfTheEarth = 0;
  else if (a == 6371229) shapeOfTheEarth = 6;
  else if (a == 6371200) shapeOfTheEarth = 8;
  else if (a == 6378160 && b == 6356775 && rf == 297) shapeOfTheEarth = 2;
  else if (a == 6378137 && b == 6356752 && rf == 298) shapeOfTheEarth = 4;
  // clang-format on

  return shapeOfTheEarth;
}

static void
gribapiDefGridLCC(grib_handle *gh, int editionNumber, int gridID, int uvRelativeToGrid)
{
  long xsize = (long) gridInqXsize(gridID);
  long ysize = (long) gridInqYsize(gridID);

  struct CDI_GridProjParams gpp;
  gridInqParamsLCC(gridID, &gpp);
  if (IS_EQUAL(gpp.x_0, gpp.mv) && IS_EQUAL(gpp.y_0, gpp.mv) && (IS_EQUAL(gpp.xval_0, gpp.mv) || IS_EQUAL(gpp.yval_0, gpp.mv)))
    {
      gpp.x_0 = gridInqXval(gridID, 0);
      gpp.y_0 = gridInqYval(gridID, 0);
    }
  gridVerifyProjParamsLCC(&gpp);
  if (gpp.xval_0 < 0.0) gpp.xval_0 += 360.0;
  if (gpp.lon_0 < 0.0) gpp.lon_0 += 360.0;

  bool isSouthPole = (gpp.lat_1 < 0.0);
  if (isSouthPole)
    {
      gpp.lat_1 = -gpp.lat_2;
      gpp.lat_2 = -gpp.lat_2;
    }
  int projflag = 0;
  if (isSouthPole) gribbyte_set_bit(&projflag, 1);

  double xinc = gridInqXinc(gridID);
  double yinc = gridInqYinc(gridID);
  if (IS_EQUAL(xinc, 0.0)) xinc = gridInqXincInMeter(gridID);
  if (IS_EQUAL(yinc, 0.0)) yinc = gridInqYincInMeter(gridID);

  static const char mesg[] = "lambert";
  size_t len = sizeof(mesg) - 1;
  GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

  GRIB_CHECK(my_grib_set_long(gh, "Nx", xsize), 0);
  GRIB_CHECK(my_grib_set_long(gh, "Ny", ysize), 0);
  GRIB_CHECK(my_grib_set_double(gh, "DxInMetres", fabs(xinc)), 0);
  GRIB_CHECK(my_grib_set_double(gh, "DyInMetres", fabs(yinc)), 0);
  GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", gpp.xval_0), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", gpp.yval_0), 0);
  if (editionNumber > 1) GRIB_CHECK(my_grib_set_double(gh, "LaDInDegrees", gpp.lat_1), 0);
  GRIB_CHECK(my_grib_set_double(gh, "LoVInDegrees", gpp.lon_0), 0);
  GRIB_CHECK(my_grib_set_double(gh, "Latin1InDegrees", gpp.lat_1), 0);
  GRIB_CHECK(my_grib_set_double(gh, "Latin2InDegrees", gpp.lat_2), 0);
  GRIB_CHECK(my_grib_set_long(gh, "projectionCentreFlag", projflag), 0);

  if (gpp.x_SP >= -180 && gpp.x_SP <= 360) GRIB_CHECK(my_grib_set_double(gh, "longitudeOfSouthernPoleInDegrees", gpp.x_SP), 0);
  if (gpp.y_SP >= -90 && gpp.y_SP <= 90) GRIB_CHECK(my_grib_set_double(gh, "latitudeOfSouthernPoleInDegrees", gpp.y_SP), 0);

  long shapeOfTheEarth = encode_shapeOfTheEarth(&gpp);
  if (shapeOfTheEarth) GRIB_CHECK(my_grib_set_long(gh, "shapeOfTheEarth", shapeOfTheEarth), 0);
  if (shapeOfTheEarth == 1) GRIB_CHECK(my_grib_set_long(gh, "radiusOfTheEarth", gpp.a), 0);

  long earthIsOblate = (shapeOfTheEarth == 2 || shapeOfTheEarth == 3 || shapeOfTheEarth == 4);
  if (earthIsOblate) GRIB_CHECK(my_grib_set_long(gh, "earthIsOblate", earthIsOblate), 0);

  if (uvRelativeToGrid >= 0) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);

  GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", (xinc < 0)), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", (yinc > 0)), 0);
}

static void
gribapiDefGridSTERE(grib_handle *gh, int gridID, int uvRelativeToGrid)
{
  long xsize = (long) gridInqXsize(gridID);
  long ysize = (long) gridInqYsize(gridID);

  struct CDI_GridProjParams gpp;
  gridInqParamsSTERE(gridID, &gpp);
  gridVerifyProjParamsSTERE(&gpp);
  if (gpp.xval_0 < 0.0) gpp.xval_0 += 360.0;
  int projflag = 0;

  double xinc = gridInqXinc(gridID);
  double yinc = gridInqYinc(gridID);

  static const char mesg[] = "polar_stereographic";
  size_t len = sizeof(mesg) - 1;
  GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

  GRIB_CHECK(my_grib_set_long(gh, "Nx", xsize), 0);
  GRIB_CHECK(my_grib_set_long(gh, "Ny", ysize), 0);
  GRIB_CHECK(my_grib_set_double(gh, "DxInMetres", xinc), 0);
  GRIB_CHECK(my_grib_set_double(gh, "DyInMetres", yinc), 0);
  GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", gpp.xval_0), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", gpp.yval_0), 0);
  GRIB_CHECK(my_grib_set_double(gh, "LaDInDegrees", gpp.lat_1), 0);
  GRIB_CHECK(my_grib_set_double(gh, "orientationOfTheGridInDegrees", gpp.lon_0), 0);
  long southPoleOnProjectionPlane = IS_EQUAL(gpp.lat_0, -90.0);
  GRIB_CHECK(my_grib_set_double(gh, "southPoleOnProjectionPlane", southPoleOnProjectionPlane), 0);
  GRIB_CHECK(my_grib_set_long(gh, "projectionCentreFlag", projflag), 0);

  long shapeOfTheEarth = encode_shapeOfTheEarth(&gpp);
  if (shapeOfTheEarth) GRIB_CHECK(my_grib_set_long(gh, "shapeOfTheEarth", shapeOfTheEarth), 0);
  if (shapeOfTheEarth == 1) GRIB_CHECK(my_grib_set_long(gh, "radiusOfTheEarth", gpp.a), 0);

  long earthIsOblate = (shapeOfTheEarth == 2 || shapeOfTheEarth == 3 || shapeOfTheEarth == 4);
  if (earthIsOblate) GRIB_CHECK(my_grib_set_long(gh, "earthIsOblate", earthIsOblate), 0);

  if (uvRelativeToGrid >= 0) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);

  GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", (xinc < 0)), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", (yinc > 0)), 0);
}

static void
gribapiDefGridHEALPIX(grib_handle *gh, int gridID, int uvRelativeToGrid)
{
  struct CDI_GridProjParams gpp;
  gridInqParamsHEALPIX(gridID, &gpp);
  gridVerifyProjParamsHEALPIX(&gpp);
  // if (gpp.xval_0 < 0.0) gpp.xval_0 += 360.0;

  static const char mesg[] = "healpix";
  size_t len = sizeof(mesg) - 1;
  GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);

  GRIB_CHECK(my_grib_set_long(gh, "Nside", gpp.nside), 0);
  GRIB_CHECK(my_grib_set_long(gh, "ordering", gpp.order), 0);
  GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", 45.0), 0);
  // GRIB_CHECK(my_grib_set_double(gh, "longitudeOfFirstGridPointInDegrees", gpp.xval_0), 0);
  /*
  long shapeOfTheEarth = encode_shapeOfTheEarth(&gpp);
  if (shapeOfTheEarth) GRIB_CHECK(my_grib_set_long(gh, "shapeOfTheEarth", shapeOfTheEarth), 0);
  if (shapeOfTheEarth == 1) GRIB_CHECK(my_grib_set_long(gh, "radiusOfTheEarth", gpp.a), 0);

  long earthIsOblate = (shapeOfTheEarth == 2 || shapeOfTheEarth == 3 || shapeOfTheEarth == 4);
  if (earthIsOblate) GRIB_CHECK(my_grib_set_long(gh, "earthIsOblate", earthIsOblate), 0);
  */
  if (uvRelativeToGrid >= 0) GRIB_CHECK(my_grib_set_long(gh, "uvRelativeToGrid", uvRelativeToGrid), 0);
}

static void
gribapiDefGridGME(grib_handle *gh, int gridID, long gridsize)
{
  GRIB_CHECK(my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_GME), 0);

  int nd = 0, ni = 0, ni2 = 0, ni3 = 0;
  gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
  GRIB_CHECK(my_grib_set_long(gh, "nd", nd), 0);
  GRIB_CHECK(my_grib_set_long(gh, "Ni", ni), 0);
  GRIB_CHECK(my_grib_set_long(gh, "n2", ni2), 0);
  GRIB_CHECK(my_grib_set_long(gh, "n3", ni3), 0);
  GRIB_CHECK(my_grib_set_long(gh, "latitudeOfThePolePoint", 90000000), 0);
  GRIB_CHECK(my_grib_set_long(gh, "longitudeOfThePolePoint", 0), 0);

  GRIB_CHECK(my_grib_set_long(gh, "numberOfDataPoints", gridsize), 0);
  GRIB_CHECK(my_grib_set_long(gh, "totalNumberOfGridPoints", gridsize), 0);
}

static void
gribapiDefGridUnstructured(grib_handle *gh, int gridID)
{
  static bool warning = true;

  int status = my_grib_set_long(gh, "gridDefinitionTemplateNumber", GRIB2_GTYPE_UNSTRUCTURED);
  if (status != 0 && warning)
    {
      warning = false;
      Warning("Can't write reference grid!");
      Warning("gridDefinitionTemplateNumber %d not found (grib2/template.3.%d.def)!", GRIB2_GTYPE_UNSTRUCTURED,
              GRIB2_GTYPE_UNSTRUCTURED);
    }
  else
    {
      int errCount = 0;
      int number = 0;
      status = cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, &number);
      if (status < 0) errCount++;
      if (number < 0) number = 0;
      GRIB_CHECK(my_grib_set_long(gh, "numberOfGridUsed", number), 0);

      int position = 0;
      status = cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, &position);
      if (status < 0) errCount++;
      if (position < 0) position = 0;
      GRIB_CHECK(my_grib_set_long(gh, "numberOfGridInReference", position), 0);

      unsigned char uuid[CDI_UUID_SIZE];
      size_t len = CDI_UUID_SIZE;
      memset(uuid, 0, len);
      int length = CDI_UUID_SIZE;
      status = cdiInqKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
      if (status < 0) errCount++;
      if (grib_set_bytes(gh, "uuidOfHGrid", uuid, &len) != 0) Warning("Can't write UUID!");

      if (warning && errCount > 0)
        {
          warning = false;
          char uuidStr[uuidNumHexChars + 1] = { 0 };
          cdiUUID2Str(uuid, uuidStr);
          Warning("GRIB2 grid parameter missing: numberOfGridUsed=%d numberOfGridInReference=%d uuidOfHGrid=%s", number, position,
                  uuidStr);
        }
    }
}

static void
gribapiDefGridSpectral(grib_handle *gh, int gridID)
{
  int trunc = gridInqTrunc(gridID);
  enum
  {
    numTruncAtt = 3
  };
  static const char truncAttNames[numTruncAtt][2] = { "J", "K", "M" };
  for (size_t i = 0; i < numTruncAtt; ++i) GRIB_CHECK(my_grib_set_long(gh, truncAttNames[i], trunc), 0);

  if (gridInqComplexPacking(gridID))
    {
      static const char truncAttNames2[numTruncAtt][3] = { "JS", "KS", "MS" };
      for (size_t i = 0; i < numTruncAtt; ++i) GRIB_CHECK(my_grib_set_long(gh, truncAttNames2[i], 20), 0);
    }
}

static void
gribapiDefPackingType(grib_handle *gh, bool lieee, bool lspectral, bool lcomplex, int comptype, size_t gridsize)
{
  static const char mesg_spectral_complex[] = "spectral_complex";
  static const char mesg_spectral_simple[] = "spectral_simple";
  static const char mesg_grid_jpeg[] = "grid_jpeg";
  static const char mesg_grid_ccsds[] = "grid_ccsds";
  static const char mesg_ieee[] = "grid_ieee";
  static const char mesg_simple[] = "grid_simple";
  const char *mesg = mesg_simple;

  if (lspectral)
    {
      mesg = lcomplex ? mesg_spectral_complex : mesg_spectral_simple;
    }
  else if (comptype == CDI_COMPRESS_JPEG && gridsize > 1)
    {
      mesg = mesg_grid_jpeg;
    }
  else if ((comptype == CDI_COMPRESS_SZIP || comptype == CDI_COMPRESS_AEC) && gridsize > 1)
    {
      mesg = mesg_grid_ccsds;
    }
  else if (lieee)
    {
      mesg = mesg_ieee;
    }

  size_t len = strlen(mesg);
  GRIB_CHECK(my_grib_set_string(gh, "packingType", mesg, &len), 0);
}

static void
gribapiDefGrid(int editionNumber, grib_handle *gh, int gridID, int comptype, int datatype, int uvRelativeToGrid)
{
  size_t gridsize = gridInqSize(gridID);
  bool gridIsRotated = false;
  bool gridIsCurvilinear = false;
  int gridtype = grbGetGridtype(&gridID, gridsize, &gridIsRotated, &gridIsCurvilinear);

  bool lieee = (editionNumber == 2 && (datatype == CDI_DATATYPE_FLT32 || datatype == CDI_DATATYPE_FLT64));
  bool lspectral = (gridtype == GRID_SPECTRAL);
  bool lcomplex = (lspectral && gridInqComplexPacking(gridID));

  if (lieee) comptype = 0;
  if (lspectral) lieee = false;

  if (lspectral)  // gridType needs to be defined before packingType !!!
    {
      static const char mesg[] = "sh";
      size_t len = sizeof(mesg) - 1;
      GRIB_CHECK(my_grib_set_string(gh, "gridType", mesg, &len), 0);
    }

  gribapiDefPackingType(gh, lieee, lspectral, lcomplex, comptype, gridsize);

  if (lieee) GRIB_CHECK(my_grib_set_long(gh, "precision", datatype == CDI_DATATYPE_FLT64 ? 2 : 1), 0);

  if (editionNumber == 2) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", (long) gridsize), 0);

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
        gribapiDefGridRegular(gh, gridID, gridtype, gridIsRotated, gridIsCurvilinear, uvRelativeToGrid);
        break;
      }
    case CDI_PROJ_LCC:
      {
        gribapiDefGridLCC(gh, editionNumber, gridID, uvRelativeToGrid);
        break;
      }
    case CDI_PROJ_STERE:
      {
        gribapiDefGridSTERE(gh, gridID, uvRelativeToGrid);
        break;
      }
    case CDI_PROJ_HEALPIX:
      {
        if (editionNumber <= 1) Error("HEALPix grid can't be stored in GRIB edition %d!", editionNumber);
        gribapiDefGridHEALPIX(gh, gridID, uvRelativeToGrid);
        break;
      }
    case GRID_SPECTRAL:
      {
        gribapiDefGridSpectral(gh, gridID);
        break;
      }
    case GRID_GME:
      {
        if (editionNumber <= 1) Error("GME grid can't be stored in GRIB edition %d!", editionNumber);
        gribapiDefGridGME(gh, gridID, (long) gridsize);
        break;
      }
    case GRID_UNSTRUCTURED:
      {
        if (editionNumber <= 1) Error("Unstructured grid can't be stored in GRIB edition %d!", editionNumber);
        gribapiDefGridUnstructured(gh, gridID);
        break;
      }
    default:
      {
        Error("Unsupported grid type: %s", gridNamePtr(gridtype));
        break;
      }
    }
}

static void
getLevelFactor(double level, long *factor, long *out_scaled_value)
{
  const double eps = 1.0e-7;

  double scaled_value = level;
  long iscaled_value = lround(scaled_value);

  long i;
  for (i = 0; (iscaled_value < (4294967295 / 10)) && (fabs(scaled_value - (double) iscaled_value) >= eps) && i < 7; i++)
    {
      scaled_value *= 10.0;
      iscaled_value = lround(scaled_value);
    }

  (*factor) = i;
  (*out_scaled_value) = iscaled_value;
}

static void
gribapiDefLevelType(grib_handle *gh, int gcinit, const char *keyname, long leveltype)
{
  bool lset = false;
  if ((leveltype == GRIB1_LTYPE_ISOBARIC_PA || leveltype == 99 || leveltype == 100) && gribEditionNumber(gh) == 1)
    {
      if (gribGetLong(gh, "indicatorOfTypeOfLevel") != leveltype) lset = true;
    }

  if (!gcinit || lset) GRIB_CHECK(my_grib_set_long(gh, keyname, leveltype), 0);
}

static void
grib1DefLevel(grib_handle *gh, int gcinit, long leveltype1, long leveltype2, bool hasBounds, double level, double dlevel1,
              double dlevel2)
{
  (void) leveltype2;
  gribapiDefLevelType(gh, gcinit, "indicatorOfTypeOfLevel", leveltype1);

  if (hasBounds)
    {
      GRIB_CHECK(my_grib_set_long(gh, "topLevel", lround(dlevel1)), 0);
      GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", lround(dlevel2)), 0);
    }
  else
    {
      GRIB_CHECK(my_grib_set_long(gh, "level", lround(level)), 0);
    }
}

static void
grib2DefLevel(grib_handle *gh, int gcinit, long leveltype1, long leveltype2, bool hasBounds, double level, double dlevel1,
              double dlevel2)
{
  gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", leveltype1);
  if (hasBounds) gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", leveltype2);

  if (!hasBounds) dlevel1 = level;

  long scaled_level, factor;
  getLevelFactor(dlevel1, &factor, &scaled_level);
  GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfFirstFixedSurface", factor), 0);
  GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfFirstFixedSurface", scaled_level), 0);

  if (hasBounds)
    {
      getLevelFactor(dlevel2, &factor, &scaled_level);
      GRIB_CHECK(my_grib_set_long(gh, "scaleFactorOfSecondFixedSurface", factor), 0);
      GRIB_CHECK(my_grib_set_long(gh, "scaledValueOfSecondFixedSurface", scaled_level), 0);
    }
}

static void
gribapiDefLevel(int editionNumber, grib_handle *gh, int zaxisID, int levelID, int gcinit, int proddef_template_num)
{
  int zaxistype = zaxisInqType(zaxisID);
  int ltype = 0, ltype2 = -1;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFSECONDFIXEDSURFACE, &ltype2);

  bool hasBounds = (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL));
  double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) : levelID + 1;
  double dlevel1 = hasBounds ? zaxisInqLbound(zaxisID, levelID) : level;
  double dlevel2 = hasBounds ? zaxisInqUbound(zaxisID, levelID) : 0.0;

  if (zaxistype == ZAXIS_GENERIC && ltype == 0)
    {
      Warning("Changed zaxis type from %s to %s", zaxisNamePtr(zaxistype), zaxisNamePtr(ZAXIS_PRESSURE));
      zaxistype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, zaxistype);
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, "Pa");
    }

  long grib_ltype;
  {
    int (*ltypeMap)(int grib_ltype) = editionNumber <= 1 ? zaxisTypeToGrib1ltype : zaxisTypeToGrib2ltype;
    grib_ltype = ltypeMap(zaxistype);
  }
  long grib_ltype2 = (ltype != ltype2 && ltype2 != -1) ? ltype2 : grib_ltype;

  void (*defLevel)(grib_handle * gh, int gcinit, long leveltype1, long leveltype2, bool hasBounds, double level, double dlevel1,
                   double dlevel2)
      = (editionNumber <= 1) ? grib1DefLevel : grib2DefLevel;

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_HEIGHT:
    case ZAXIS_ALTITUDE:
    case ZAXIS_SIGMA:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
        if (zaxistype == ZAXIS_HEIGHT)
          {
            double sf = zaxis_units_to_meter(zaxisID);
            level *= sf;
            dlevel1 *= sf;
            dlevel2 *= sf;
          }

        /* GRIB2: PRODUCT DEFINITION TEMPLATE NUMBER 32:

           "Analysis or forecast at a horizontal level or in a
           horizontal layer at a point in time for simulate
           (synthetic) satellite data"

           The key/value pairs that are set in "grib2DefLevel" do not
           exist for this template. */
        if (editionNumber <= 1 || proddef_template_num != 32)
          defLevel(gh, gcinit, grib_ltype, grib_ltype2, hasBounds, level, dlevel1, dlevel2);

        break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TROPOPAUSE:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_LAKE_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM:
    case ZAXIS_SEDIMENT_BOTTOM_TA:
    case ZAXIS_SEDIMENT_BOTTOM_TW:
    case ZAXIS_MIX_LAYER:
    case ZAXIS_ATMOSPHERE:
      {
        defLevel(gh, gcinit, grib_ltype, grib_ltype, hasBounds, level, dlevel1, dlevel2);
      }
      break;
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        if (editionNumber <= 1)
          {
            grib_ltype = hasBounds ? GRIB1_LTYPE_HYBRID_LAYER : GRIB1_LTYPE_HYBRID;
          }
        defLevel(gh, gcinit, grib_ltype, grib_ltype, hasBounds, level, dlevel1, dlevel2);

        if (!gcinit)
          {
            int vctsize = zaxisInqVctSize(zaxisID);
            if (vctsize > 0)
              {
                GRIB_CHECK(my_grib_set_long(gh, "PVPresent", 1), 0);
                GRIB_CHECK(grib_set_double_array(gh, "pv", zaxisInqVctPtr(zaxisID), (size_t) vctsize), 0);
              }
          }

        break;
      }
    case ZAXIS_PRESSURE:
      {
        if (level < 0) Warning("Pressure level of %f Pa is below zero!", level);

        if (!zaxis_units_is_Pa(zaxisID))
          {
            level *= 100;
            dlevel1 *= 100;
            dlevel2 *= 100;
          }

        if (editionNumber <= 1)
          {
            double dum;
            if (level < 32768 && (level < 100 || modf(level / 100, &dum) > 0))
              grib_ltype = GRIB1_LTYPE_ISOBARIC_PA;
            else
              level /= 100;
          }
        else if (ltype2 == -1)
          ltype2 = GRIB2_LTYPE_ISOBARIC;
        defLevel(gh, gcinit, grib_ltype, ltype2, hasBounds, level, dlevel1, dlevel2);

        break;
      }
    case ZAXIS_SNOW:
      if (editionNumber <= 1)
        ;  // not available
      else
        {
          grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype, hasBounds, level, dlevel1, dlevel2);
        }

      break;
    case ZAXIS_DEPTH_BELOW_LAND:
      {
        double sf = editionNumber <= 1 ? zaxis_units_to_centimeter(zaxisID) : zaxis_units_to_meter(zaxisID);
        grib_ltype = editionNumber <= 1 ? (hasBounds ? GRIB1_LTYPE_LANDDEPTH_LAYER : GRIB1_LTYPE_LANDDEPTH) : grib_ltype;
        defLevel(gh, gcinit, grib_ltype, grib_ltype, hasBounds, level * sf, dlevel1 * sf, dlevel2 * sf);

        break;
      }
    case ZAXIS_REFERENCE:
      {
        if (!gcinit) GRIB_CHECK(my_grib_set_long(gh, "genVertHeightCoords", 1), 0);

        if (editionNumber <= 1)
          ;  // not available
        else
          {
            if (hasBounds)
              {
                gribapiDefLevelType(gh, gcinit, "typeOfFirstFixedSurface", grib_ltype);
                gribapiDefLevelType(gh, gcinit, "typeOfSecondFixedSurface", grib_ltype2);
                GRIB_CHECK(my_grib_set_long(gh, "topLevel", (long) dlevel1), 0);
                GRIB_CHECK(my_grib_set_long(gh, "bottomLevel", (long) dlevel2), 0);
              }
            else
              {
                grib2DefLevel(gh, gcinit, grib_ltype, grib_ltype2, hasBounds, level, dlevel1, dlevel2);
              }

            GRIB_CHECK(my_grib_set_long(gh, "NV", 6), 0);
            int number = 0;
            cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, &number);
            GRIB_CHECK(my_grib_set_long(gh, "numberOfVGridUsed", number), 0);
            int nlev = 0;
            cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NLEV, &nlev);
            GRIB_CHECK(my_grib_set_long(gh, "nlev", nlev), 0);
            unsigned char uuid[CDI_UUID_SIZE];
            int length = CDI_UUID_SIZE;
            memset(uuid, 0, length);
            cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuid, &length);
            size_t len = CDI_UUID_SIZE;
            if (grib_set_bytes(gh, "uuidOfVGrid", uuid, &len) != 0) Warning("Can't write UUID!");
          }

        break;
      }
    case ZAXIS_GENERIC:
      {
        defLevel(gh, gcinit, ltype, ltype, hasBounds, level, dlevel1, dlevel2);
        break;
      }
    default:
      {
        Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
        break;
      }
    }
}

int
gribapiGetScanningMode(grib_handle *gh)
{
  long iScansNegatively, jScansPositively, jPointsAreConsecutive;
  GRIB_CHECK(grib_get_long(gh, "iScansNegatively", &iScansNegatively), 0);
  GRIB_CHECK(grib_get_long(gh, "jScansPositively", &jScansPositively), 0);
  GRIB_CHECK(grib_get_long(gh, "jPointsAreConsecutive", &jPointsAreConsecutive), 0);
  int scanningMode = 128 * (bool) iScansNegatively + 64 * (bool) jScansPositively + 32 * (bool) jPointsAreConsecutive;
  if (cdiDebugExt >= 30)
    printf("gribapiGetScanningMode(): Scanning mode = %02d (%1d%1d%1d)*32; \n", scanningMode, (int) jPointsAreConsecutive,
           (int) jScansPositively, (int) iScansNegatively);

  return scanningMode;
}

void
gribapiSetScanningMode(grib_handle *gh, int scanningMode)
{
  // 127: reserved for testing; generated test data will be in 64 scanning mode
  // if (scanningMode== 127)  scanningMode = 64;

  long iScansNegatively = (scanningMode & 128) / 128;
  long jScansPositively = (scanningMode & 64) / 64;
  long jPointsAreConsecutive = (scanningMode & 32) / 32;

  if (cdiDebugExt >= 30 && gribEditionNumber(gh) <= 1)
    {
      long paramId, levelTypeId, levelId, uvRelativeToGrid;
      GRIB_CHECK(grib_get_long(gh, "uvRelativeToGrid", &uvRelativeToGrid), 0);
      GRIB_CHECK(grib_get_long(gh, "indicatorOfParameter", &paramId), 0);
      GRIB_CHECK(grib_get_long(gh, "indicatorOfTypeOfLevel", &levelTypeId), 0);
      GRIB_CHECK(grib_get_long(gh, "level", &levelId), 0);
      printf("gribapiSetScanningMode(): (param,ltype,level) = (%3d,%3d,%4d); Scanning mode = %02d (%1d%1d%1d)*32;  "
             "uvRelativeToGrid = %02d\n",
             (int) paramId, (int) levelTypeId, (int) levelId, scanningMode, (int) jPointsAreConsecutive, (int) jScansPositively,
             (int) iScansNegatively, (int) uvRelativeToGrid);
    }

  GRIB_CHECK(my_grib_set_long(gh, "iScansNegatively", iScansNegatively), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jScansPositively", jScansPositively), 0);
  GRIB_CHECK(my_grib_set_long(gh, "jPointsAreConsecutive", jPointsAreConsecutive), 0);
}

/*
  TABLE 8. SCANNING MODE FLAG

  (GDS Octet 28)
  BIT     VALUE     MEANING
  1       0       Points scan in +i direction
          1       Points scan in -i direction
  2       0       Points scan in -j direction
          1       Points scan in +j direction
  3       0       Adjacent points in i direction are consecutive
                    (FORTRAN: (I,J))
          1       Adjacent points in j direction are consecutive
                  (FORTRAN: (J,I))

  => Scanning Mode     0 0 0 0 0 0 0 0  (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & North->South)
  => Scanning Mode     0 1 0 0 0 0 0 0  (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & South->North )
  => Scanning Mode     1 1 0 0 0 0 0 0  (96 dec)  +i, +j; j direction consecutive (column-major order South->North & West->East )

  NOTE:  South->North  - As if you would plot the data as image on the screen
                         where [0,0] of the data is the top-left pixel.

                         grib2ppm LAMH_D11_201302150000_00000_oro | display ppm:-
                         ImageMagick (display): [0,0] of an image belongs to the top-left pixel
  [DEFAULT] : 64 dec

  iScansNegatively = 0;
  jScansPositively = 1;
  jPointsAreConsecutive = 0;    => Scanning Mode 64

  cdo selindexbox,1,726,100,550 LAMH_D11_201302150000_00000_oro LAMH_D11_201302150000_00000_oro_cropped
  grib2ppm LAMH_D11_201302150000_00000_oro_cropped | /usr/bin/display ppm:- &
  # ^^^ this image will be missing the souther parts of data

  grib2ppm LAMH_D11_201302150000_00000_oro | /usr/bin/display ppm:- &
  # ^ full domain data
*/

#ifdef HIRLAM_EXTENSIONS
static void
verticallyFlipGridDefinitionWhenScanningModeChanged(grib_handle *gh, double yfirst, double ylast, double yinc)
{
  /*
  Nj = 550;
  latitudeOfFirstGridPointInDegrees = -30.8;
  latitudeOfLastGridPointInDegrees = 24.1;
  iScansNegatively = 0;
  jScansPositively = 0;
  jPointsAreConsecutive = 0;
  jDirectionIncrementInDegrees = 0.1;

  When switching from scanning mode 0 <=> 64
  yfirst = -30.8 + (550-1)*0.1

  yfirst = yfirst + (ysize-1) * yinc
  yinc   = -1.0*yinc
  */

  // long jDim=0;
  // GRIB_CHECK(grib_get_long(gh, "Nj", &jDim), 0);

  double latitudeOfFirstGridPointInDegrees;
  double latitudeOfLastGridPointInDegrees;
  double jDirectionIncrementInDegrees;

  // GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &latitudeOfFirstGridPointInDegrees), 0);  // yfirst
  // GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees", &latitudeOfLastGridPointInDegrees), 0);    // ylast
  // GRIB_CHECK(grib_get_double(gh, "jDirectionIncrementInDegrees", &jDirectionIncrementInDegrees), 0);  // yinc

  if (cdiDebugExt >= 10) Message(" BEFORE: yfirst = %f; ylast = %f; yinc = %f; ", yfirst, ylast, yinc);

  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfFirstGridPointInDegrees", ylast), 0);
  GRIB_CHECK(my_grib_set_double(gh, "latitudeOfLastGridPointInDegrees", yfirst), 0);
  // yinc *= -1.0; // don't set yinc here ...
  // GRIB_CHECK(my_grib_set_double(gh, "jDirectionIncrementInDegrees", yinc), 0);

  if (cdiDebugExt >= 10)
    {
      GRIB_CHECK(grib_get_double(gh, "latitudeOfFirstGridPointInDegrees", &latitudeOfFirstGridPointInDegrees), 0);  // yfirst
      GRIB_CHECK(grib_get_double(gh, "latitudeOfLastGridPointInDegrees", &latitudeOfLastGridPointInDegrees), 0);    // ylast
      GRIB_CHECK(grib_get_double(gh, "jDirectionIncrementInDegrees", &jDirectionIncrementInDegrees), 0);            // yinc
      Message("CHANGED INTO:  yfirst = %f, ylast = %f, yinc = %f", latitudeOfFirstGridPointInDegrees,
              latitudeOfLastGridPointInDegrees, jDirectionIncrementInDegrees);
    }
}

static void
convertDataScanningMode(int scanModeIN, int scanModeOUT, double *data, size_t gridsize, size_t iDim, size_t jDim)
{
  size_t idxIN, idxOUT;

  // 127: reserved for testing; it will generate test data in 64 scanning mode
  if (scanModeOUT == 127)  // fill with testdata ...
    {
      scanModeOUT = 64;
      if (cdiDebugExt >= 30) printf("convertDataScanningMode(): Generating test data in 64 scanning mode..\n");
      for (size_t j = 0; j < jDim; j++)
        {
          size_t jXiDim = j * iDim;
          for (size_t i = 0; i < iDim; i++)
            {
              idxIN = i + jXiDim;
              data[idxIN] = (double) (100.0 * j + i);
            }
        }
    }

  if ((iDim * jDim) != gridsize)
    {
      if (cdiDebugExt >= 30)
        printf("convertDataScanningMode(): ERROR: (iDim*jDim)!= gridsize;  (%zu * %zu) != %zu\n", iDim, jDim, gridsize);
      return;
    }
  if (cdiDebugExt >= 30)
    printf("convertDataScanningMode(): scanModeIN=%02d => scanModeOUT=%02d ; where: (iDim * jDim == gridsize)  (%zu*%zu == %zu)\n",
           scanModeIN, scanModeOUT, iDim, jDim, gridsize);

  if (cdiDebugExt >= 100)
    {
      printf("convertDataScanningMode(): data IN:\n");
      for (size_t j = 0; j < jDim; j++)
        {
          size_t jXiDim = j * iDim;
          for (size_t i = 0; i < iDim; i++)
            {
              idxIN = i + jXiDim;
              printf("%03.0f, ", data[idxIN]);
            }
          printf("\n");
        }
    }

  if (scanModeIN == scanModeOUT)
    {
      if (cdiDebugExt >= 30) printf("convertDataScanningMode(): INFO: Nothing to do;  scanModeIN==scanModeOUT..\n");
      return;
    }

  if (0)
    {
      return;
      if (scanModeOUT == 00)
        {
          if (cdiDebugExt > 0) printf("convertDataScanningMode(): Leave data unchaged BUT set scanModeOUT=00.\n");
          // CHECK:  Looks like that GRIB-API provide (no matter what) data in the scannning mode 00, even it is store in the
          // gribfile as 64 !!
          return;
        }
    }
  double *dataCopy = (double *) Malloc(gridsize * sizeof(double));
  memcpy((void *) dataCopy, (void *) data, gridsize * sizeof(double));

  if (scanModeIN
      == 64)  // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
    {         // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
              // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
      if (scanModeOUT == 00)
      // CHECK:  Looks like that GRIB-API provide (no matter what) data in the scannning mode 00, even it is store in the gribfile
      // as 64 !!
#define VERTICAL_FLIP
#ifdef VERTICAL_FLIP
        {  // flip the data vertically ..
          idxIN = 0;
          idxOUT = (jDim - 1) * iDim;
          if (cdiDebugExt >= 30) printf("convertDataScanningMode():  copying rows nr. (%04d : %04zu)\n", 0, jDim - 1);
          for (size_t j = 0; j < jDim; j++)
            {
              memcpy((void *) &data[idxOUT], (void *) &dataCopy[idxIN], iDim * sizeof(double));
              idxIN += iDim;
              idxOUT -= iDim;
            }
        }  // end if (scanModeOUT==00)*/
#endif
#ifdef HORIZONTAL_FLIP
      {  // flip data horizontally ...
        if (1)
          {
            if (cdiDebugExt >= 30) printf("convertDataScanningMode():  copying columns nr. (%04d : %04d);\n", 0, iDim - 1);
            for (size_t i = 0; i < iDim; i++)
              {
                for (size_t j = 0; j < jDim; j++)
                  {
                    size_t jXiDim = j * iDim;
                    idxIN = i + jXiDim;
                    // data[idxIN] = (double) (100.0*j +i);  // just some testdata ..
                    idxOUT = iDim - i - 1 + jXiDim;
                    // printf("[%03d=>%03d] = %f;",idxIN,idxOUT,dataCopy[idxIN]);
                    data[idxOUT] = dataCopy[idxIN];
                  }
              }
          }
      }  // end if (scanModeOUT==00)
#endif

      if (scanModeOUT == 96)
        {  // transpose the data
          if (cdiDebugExt >= 30)
            printf("convertDataScanningMode():  transpose data rows=>columns nr. (%04d : %04zu) => (%04d : %04zu);\n", 0, iDim - 1,
                   0, jDim - 1);
          for (size_t j = 0; j < jDim; j++)
            {
              size_t jXiDim = j * iDim;
              for (size_t i = 0; i < iDim; i++)
                {
                  idxIN = i + jXiDim;
                  idxOUT = j + i * jDim;
                  // printf("[%03d=>%03d] = %f;",idxIN,idxOUT,dataCopy[idxIN]);
                  data[idxOUT] = dataCopy[idxIN];
                }
              // printf(".\n");
            }
        }  // end if (scanModeOUT==96)
    }      // end if (scanModeIN==64)

  if (scanModeIN
      == 00)  // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
    {         // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
              // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
      if (scanModeOUT == 64)
        {  // flip the data vertically ..
          idxIN = 0;
          idxOUT = (jDim - 1) * iDim;
          for (size_t j = 0; j < jDim; j++)
            {
              if (cdiDebugExt >= 25)
                printf("convertDataScanningMode():  copying row nr. %04zu; [idxIN=%08zu] => [idxOUT=%08zu]\n", j, idxIN, idxOUT);
              memcpy((void *) &data[idxOUT], (void *) &dataCopy[idxIN], iDim * sizeof(double));
              idxIN += iDim;
              idxOUT -= iDim;
            }
        }  // end if (scanModeOUT==64)

      if (scanModeOUT == 96)
        {  // transpose the data
          size_t jInv;
          for (size_t j = 0; j < jDim; j++)
            {
              if (cdiDebugExt >= 30) printf("convertDataScanningMode():  processing row nr. %04zu;\n", j);
              jInv = (jDim - 1) - j;
              for (size_t i = 0; i < iDim; i++) data[j + i * jDim] = dataCopy[i + jInv * iDim];  // source data has -j
            }
        }  // end if (scanModeOUT==96)
    }      // end if (scanModeIN==00)

  if (scanModeIN
      == 96)  // Scanning Mode (00 dec)  +i, -j; i direction consecutive (row-major    order West->East   & South->North )
    {         // Scanning Mode (64 dec)  +i, +j; i direction consecutive (row-major    order West->East   & North->South )
              // Scanning Mode (96 dec)  +i, +j; j direction consecutive (column-major order North->South & West->East )
      if (scanModeOUT == 64)
        {  // transpose the data
          for (size_t j = 0; j < jDim; j++)
            {
              if (cdiDebugExt >= 30) printf("convertDataScanningMode():  processing row nr. %04zu;\n", j);
              size_t jXiDim = j * iDim;
              for (size_t i = 0; i < iDim; i++)
                // data[j + i*jDim] =  dataCopy[i + j*iDim];
                data[i + jXiDim] = dataCopy[j + i * jDim];
            }
        }  // end if (scanModeOUT==64)

      if (scanModeOUT == 00)
        {  // transpose the data
          idxIN = 0;
          idxOUT = 0;
          size_t jInv;
          for (size_t j = 0; j < jDim; j++)
            {
              if (cdiDebugExt >= 30) printf("convertDataScanningMode():  processing row nr. %04zu;\n", j);
              jInv = (jDim - 1) - j;
              size_t jXiDim = j * iDim;
              for (size_t i = 0; i < iDim; i++)
                // data[jInv + iXjDim] =  dataCopy[i + jXiDim];  // target data has -j
                data[i + jXiDim] = dataCopy[jInv + i * jDim];  // target data has -j
            }
        }  // end if (scanModeOUT==00)
    }      // end if (scanModeIN==96)

  if (cdiDebugExt >= 100)
    {
      printf("convertDataScanningMode(): data OUT (new scanning mode):\n");
      for (size_t j = 0; j < jDim; j++)
        {
          size_t jXiDim = j * iDim;
          for (size_t i = 0; i < iDim; i++)
            {
              idxIN = i + jXiDim;
              printf("%03.0f, ", data[idxIN]);
            }
          printf("\n");
        }
    }

  free(dataCopy);
}
#endif  // HIRLAM_EXTENSIONS

static void
gribapiSetExtMode(grib_handle *gh, int gridID, size_t datasize, const void *data)
{
#ifndef HIRLAM_EXTENSIONS
  (void) gh;
  (void) data;
  (void) datasize;
#endif
  int gridtype = gridInqType(gridID);
  if (gridtype == GRID_PROJECTION)
    {
      int projtype = gridInqProjType(gridID);
      // clang-format off
      if      (projtype == CDI_PROJ_RLL)   gridtype = GRID_LONLAT;
      else if (projtype == CDI_PROJ_LCC)   gridtype = CDI_PROJ_LCC;
      else if (projtype == CDI_PROJ_STERE) gridtype = CDI_PROJ_STERE;
      // clang-format on
    }

  if (gridtype == GRID_GENERIC || gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED
      || gridtype == CDI_PROJ_LCC)
    {
#ifdef HIRLAM_EXTENSIONS
      int scanModeIN = 0;
      cdiInqKeyInt(gridID, CDI_GLOBAL, CDI_KEY_SCANNINGMODE, &scanModeIN);

      if (cdiDebugExt >= 100) Message("scanModeIN=%d; gridsize=%zu", scanModeIN, gridInqSize(gridID));

      if (cdiGribDataScanningMode.active)  // allowed modes: <0, 64, 96>; Default is 64
        {
          size_t iDim = gridInqXsize(gridID);
          size_t jDim = gridInqYsize(gridID);

          double yfirst = gridInqYval(gridID, 0);
          double ylast = gridInqYval(gridID, jDim - 1);
          double yinc = gridInqYinc(gridID);

          int scanModeOUT = cdiGribDataScanningMode.value;
          convertDataScanningMode(scanModeIN, scanModeOUT, (double *) data, datasize, iDim, jDim);
          // This will overrule the old scanning mode of the given grid
          if (cdiDebugExt >= 10) Message("Set GribDataScanningMode (%d) => (%d)", scanModeIN, cdiGribDataScanningMode.value);
          gribapiSetScanningMode(gh, cdiGribDataScanningMode.value);

          if (((scanModeIN == 00) && (cdiGribDataScanningMode.value == 64))
              || ((scanModeIN == 64) && (cdiGribDataScanningMode.value == 00)))
            verticallyFlipGridDefinitionWhenScanningModeChanged(gh, yfirst, ylast, yinc);
        }
      else
        {
          if (cdiDebugExt >= 100) Message("Set GribDataScanningMode => (%d) based on used grid", scanModeIN);
          gribapiSetScanningMode(gh, scanModeIN);
        }
#endif
    }
}

// #define GRIBAPIENCODETEST 1

size_t
gribapiEncode(int memType, int varID, int levelID, int vlistID, int gridID, int zaxisID, CdiDateTime vDateTime, int tsteptype,
              int numavg, size_t datasize, const void *data, size_t numMissVals, void **gribbuffer, size_t *gribbuffersize,
              int comptype, void *gribContainer)
{
  long editionNumber = 2;

  // extern unsigned char _grib_template_GRIB2[];
  cdi_check_gridsize_int_limit("GRIB", datasize);

  int param = vlistInqVarParam(vlistID, varID);
  int datatype = vlistInqVarDatatype(vlistID, varID);
  int typeOfGeneratingProcess = 0;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, &typeOfGeneratingProcess);
  int productDefinitionTemplate = -1;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, &productDefinitionTemplate);

  int uvRelativeToGrid = -1;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, &uvRelativeToGrid);

#ifdef GRIBAPIENCODETEST
  grib_handle *gh = (grib_handle *) gribHandleNew(editionNumber);
#else
  gribContainer_t *gc = (gribContainer_t *) gribContainer;
  assert(gc != NULL);
  grib_handle *gh = (struct grib_handle *) gc->gribHandle;
#endif

  GRIB_CHECK(grib_get_long(gh, "editionNumber", &editionNumber), 0);

  if (editionNumber == 2)
    {
      if (!gc->init)
        {
          int backgroundProcess = 0;
          cdiInqKeyInt(vlistID, varID, CDI_KEY_BACKGROUNDPROCESS, &backgroundProcess);
          GRIB_CHECK(my_grib_set_long(gh, "typeOfGeneratingProcess", typeOfGeneratingProcess), 0);
          GRIB_CHECK(my_grib_set_long(gh, "backgroundProcess", backgroundProcess), 0);
          int status, tablesVersion, localTablesVersion;
          status = cdiInqKeyInt(vlistID, varID, CDI_KEY_TABLESVERSION, &tablesVersion);
          if (status == 0) GRIB_CHECK(my_grib_set_long(gh, "tablesVersion", (long) tablesVersion), 0);
          status = cdiInqKeyInt(vlistID, varID, CDI_KEY_LOCALTABLESVERSION, &localTablesVersion);
          if (status == 0) GRIB_CHECK(my_grib_set_long(gh, "localTablesVersion", (long) localTablesVersion), 0);
          int typeOfProcessedData = 0;
          status = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFPROCESSEDDATA, &typeOfProcessedData);
          if (status == 0) GRIB_CHECK(my_grib_set_long(gh, "typeOfProcessedData", (long) typeOfProcessedData), 0);
          /*
          int constituentType = 0;
          status = cdiInqKeyInt(vlistID, varID, CDI_KEY_CONSTITUENTTYPE, &constituentType);
          if ( status == 0 ) GRIB_CHECK(my_grib_set_long(gh, "constituentType", (long)constituentType), 0);
          */
        }
    }

  gribapiDefTime((int) editionNumber, productDefinitionTemplate, typeOfGeneratingProcess, gh, vDateTime, tsteptype, numavg,
                 vlistInqTaxis(vlistID), gc->init);

  {
    int typeOfTimeIncrement = 0;
    int status = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFTIMEINCREMENT, &typeOfTimeIncrement);
    if (status == 0) grib_set_long(gh, "typeOfTimeIncrement", (long) typeOfTimeIncrement);
  }

  {
    int status, perturbationNumber, numberOfForecastsInEnsemble, typeOfEnsembleForecast;
    status = cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble);
    if (status == 0) grib_set_long(gh, "numberOfForecastsInEnsemble", numberOfForecastsInEnsemble);
    status = cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber);
    if (status == 0) grib_set_long(gh, "perturbationNumber", perturbationNumber);
    status = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, &typeOfEnsembleForecast);
    if (status == 0) grib_set_long(gh, "typeOfEnsembleForecast", typeOfEnsembleForecast);
  }

  if (!gc->init) gribapiDefInstitut(gh, vlistID, varID);
  if (!gc->init) gribapiDefModel(gh, vlistID, varID);

  if (!gc->init)
    {
      char name[256], stdname[256];
      vlistInqVarName(vlistID, varID, name);
      vlistInqVarStdname(vlistID, varID, stdname);
      gribapiDefParam((int) editionNumber, gh, param, name, stdname);
    }

  if (!gc->init && editionNumber == 2)
    {
      int shapeOfTheEarth = 0;
      cdiInqKeyInt(vlistID, varID, CDI_KEY_SHAPEOFTHEEARTH, &shapeOfTheEarth);
      GRIB_CHECK(my_grib_set_long(gh, "shapeOfTheEarth", (long) shapeOfTheEarth), 0);

      int grib2LocalSectionNumber, section2PaddingLength;
      int mpimType, mpimClass, mpimUser;
      if (cdiInqKeyInt(vlistID, varID, CDI_KEY_MPIMTYPE, &mpimType) == CDI_NOERR
          && cdiInqKeyInt(vlistID, varID, CDI_KEY_MPIMCLASS, &mpimClass) == CDI_NOERR
          && cdiInqKeyInt(vlistID, varID, CDI_KEY_MPIMUSER, &mpimUser) == CDI_NOERR)
        {
          grib_set_long(gh, "grib2LocalSectionPresent", 1);
          grib_set_long(gh, "grib2LocalSectionNumber", 1);
          grib_set_long(gh, "mpimType", mpimType);
          grib_set_long(gh, "mpimClass", mpimClass);
          grib_set_long(gh, "mpimUser", mpimUser);

          int revNumLen = 20;
          unsigned char revNumber[revNumLen];
          if (cdiInqKeyBytes(vlistID, varID, CDI_KEY_REVNUMBER, revNumber, &revNumLen) == CDI_NOERR)
            {
              size_t revNumLenS = revNumLen;
              grib_set_bytes(gh, "revNumber", revNumber, &revNumLenS);
            }
          int revStatus;
          if (cdiInqKeyInt(vlistID, varID, CDI_KEY_REVSTATUS, &revStatus) == CDI_NOERR) grib_set_long(gh, "revStatus", revStatus);
        }
      else if (cdiInqKeyInt(vlistID, varID, CDI_KEY_GRIB2LOCALSECTIONNUMBER, &grib2LocalSectionNumber) == CDI_NOERR
               && cdiInqKeyInt(vlistID, varID, CDI_KEY_SECTION2PADDINGLENGTH, &section2PaddingLength) == CDI_NOERR)
        {
          grib_set_long(gh, "grib2LocalSectionPresent", 1);
          grib_set_long(gh, "grib2LocalSectionNumber", grib2LocalSectionNumber);
          unsigned char *section2Padding = (unsigned char *) Malloc(section2PaddingLength);
          cdiInqKeyBytes(vlistID, varID, CDI_KEY_SECTION2PADDING, section2Padding, &section2PaddingLength);
          size_t len = section2PaddingLength;
          // Does not work anymore with ecCodes 2.22.0/2.25.0!!!
          // ECCODES ERROR   :  pack_bytes: Wrong size (10) for section2Padding. It is 0 bytes long
          grib_set_bytes(gh, "section2Padding", section2Padding, &len);
          Free(section2Padding);
        }
    }

  // bitsPerValue have to be defined first (complex packing)
  GRIB_CHECK(my_grib_set_long(gh, "bitsPerValue", (long) grbBitsPerValue(datatype)), 0);

  if (!gc->init) gribapiDefGrid((int) editionNumber, gh, gridID, comptype, datatype, uvRelativeToGrid);

  gribapiDefLevel((int) editionNumber, gh, zaxisID, levelID, gc->init, productDefinitionTemplate);

  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int zaxisSize = zaxisInqSize(zaxisID);
  // if (!gc->init)
  {
    int ret = 0;

    // NOTE: Optional key/value pairs: Note that we do not distinguish between tiles here!

    for (int i = 0; i < vlistptr->vars[varID].opt_grib_nentries; i++)
      {
        if (vlistptr->vars[varID].opt_grib_kvpair[i].update)
          {
            // DR: Fix for multi-level fields (otherwise only the 1st level is correct)
            if (zaxisSize == (levelID + 1)) vlistptr->vars[varID].opt_grib_kvpair[i].update = false;

            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_double)
              {
                if (CDI_Debug)
                  Message("key \"%s\"  :   double value = %g", vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                my_grib_set_double(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                   vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val);
                GRIB_CHECK(ret, 0);
              }
            if (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_int)
              {
                if (CDI_Debug)
                  Message("key \"%s\"  :   integer value = %d", vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                          vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                my_grib_set_long(gh, vlistptr->vars[varID].opt_grib_kvpair[i].keyword,
                                 (long) vlistptr->vars[varID].opt_grib_kvpair[i].int_val);
                GRIB_CHECK(ret, 0);
              }
          }
      }
  }

  if (numMissVals > 0)
    {
      GRIB_CHECK(my_grib_set_long(gh, "bitmapPresent", 1), 0);
      GRIB_CHECK(my_grib_set_double(gh, "missingValue", vlistInqVarMissval(vlistID, varID)), 0);
    }

  gribapiSetExtMode(gh, gridID, datasize, data);

  if (memType == MEMTYPE_FLOAT)
    {
#ifdef HAVE_GRIBAPI_FLOAT_INTERFACE
      GRIB_CHECK(grib_set_float_array(gh, "values", (const float *) data, datasize), 0);
#else
      Error("grib_set_float_array() not found!");
#endif
    }
  else
    {
      GRIB_CHECK(grib_set_double_array(gh, "values", (const double *) data, datasize), 0);
    }

  if (numMissVals)
    {
      long numberOfMissing = -1;
      GRIB_CHECK(grib_get_long(gh, "numberOfMissing", &numberOfMissing), 0);
      if (numberOfMissing == 0) GRIB_CHECK(my_grib_set_long(gh, "bitmapPresent", 0), 0);
    }

  // get the size of coded message
  const void *dummy = NULL;
  size_t recsize = 0;
  GRIB_CHECK(grib_get_message(gh, &dummy, &recsize), 0);
  recsize += 512;  // add some space for possible filling
  *gribbuffersize = recsize;
  *gribbuffer = Malloc(*gribbuffersize);

  if (CDI_Debug && !gc->init && editionNumber == 2)
    {
      long pdis;
      grib_get_long(gh, "discipline", &pdis);
      if (pdis != 255)
        {
          char cdi_name[CDI_MAX_NAME];
          cdi_name[0] = 0;
          vlistInqVarName(vlistID, varID, cdi_name);
          char grb_name[256];
          gribapi_get_string(gh, "shortName", grb_name, sizeof(grb_name));
          str_to_lower(cdi_name);
          str_to_lower(grb_name);
          bool checkName = (!grb_name[0] && strncmp(cdi_name, "param", 5) == 0) ? false : true;
          if (checkName && ((strlen(cdi_name) != strlen(grb_name)) || !strStartsWith(cdi_name, grb_name)))
            Warning("*** GRIB2 shortName does not correspond to chosen variable name: \"%s\" (\"%s\").",
                    grb_name[0] ? grb_name : "unknown", cdi_name);
        }
    }

  // get a copy of the coded message
  GRIB_CHECK(grib_get_message_copy(gh, *gribbuffer, &recsize), 0);

#ifdef GRIBAPIENCODETEST
  gribHandleDelete(gh);
#endif

  gc->init = true;

  return recsize;
}

void
gribapiChangeParameterIdentification(grib_handle *gh, int code, int ltype, int level)
{
  //  timeRangeIndicator: could be included later
  if (gribEditionNumber(gh) <= 1)
    {
      if (code != -1) GRIB_CHECK(my_grib_set_long(gh, "indicatorOfParameter", code), 0);
      if (ltype != -1) GRIB_CHECK(my_grib_set_long(gh, "indicatorOfTypeOfLevel", ltype), 0);
      if (level != -1) GRIB_CHECK(my_grib_set_long(gh, "level", level), 0);
    }
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
