#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "cdf_util.h"
#include "error.h"

char *
str_to_lower(char *str)
{
  if (str)
    for (size_t i = 0; str[i]; ++i) str[i] = (char) tolower((int) str[i]);

  return str;
}

bool
strStartsWith(const char *vstr, const char *cstr)
{
  bool is_equal = false;
  if (vstr && cstr)
    {
      size_t clen = strlen(cstr);
      size_t vlen = strlen(vstr);
      if (clen <= vlen) is_equal = (memcmp(vstr, cstr, clen) == 0);
    }
  return is_equal;
}

int
get_time_units(size_t len, const char *ptu)
{
  int timeunit = -1;

  while (isspace(*ptu) && len)
    {
      ptu++;
      len--;
    }

  // clang-format off
  if (len > 2)
    {
      if      (strStartsWith(ptu, "sec"))            timeunit = TUNIT_SECOND;
      else if (strStartsWith(ptu, "minute"))         timeunit = TUNIT_MINUTE;
      else if (strStartsWith(ptu, "hour"))           timeunit = TUNIT_HOUR;
      else if (strStartsWith(ptu, "day"))            timeunit = TUNIT_DAY;
      else if (strStartsWith(ptu, "month"))          timeunit = TUNIT_MONTH;
      else if (strStartsWith(ptu, "calendar_month")) timeunit = TUNIT_MONTH;
      else if (strStartsWith(ptu, "year"))           timeunit = TUNIT_YEAR;
    }
  else if     (len == 1 && ptu[0] == 's')            timeunit = TUNIT_SECOND;
  // clang-format on

  return timeunit;
}

bool
is_time_units(const char *timeunits)
{
  while (isspace(*timeunits)) timeunits++;

  // clang-format off
  return (strStartsWith(timeunits, "sec")
       || strStartsWith(timeunits, "minute")
       || strStartsWith(timeunits, "hour")
       || strStartsWith(timeunits, "day")
       || strStartsWith(timeunits, "month")
       || strStartsWith(timeunits, "calendar_month")
       || strStartsWith(timeunits, "year"));
  // clang-format on
}

bool
is_timeaxis_units(const char *timeunits)
{
  bool status = false;

  size_t len = strlen(timeunits);
  char *tu = (char *) malloc((len + 1) * sizeof(char));

  for (size_t i = 0; i < len; i++) tu[i] = (char) tolower((int) timeunits[i]);

  int timeunit = get_time_units(len, tu);
  if (timeunit != -1)
    {
      size_t pos = 0;
      while (!isspace(tu[pos]) && tu[pos] != 0) pos++;
      if (tu[pos])
        {
          while (isspace(tu[pos])) pos++;

          status = strStartsWith(tu + pos, "as") || strStartsWith(tu + pos, "since");
        }
    }

  free(tu);

  return status;
}

bool
is_height_units(const char *units)
{
  int u0 = units[0];

  // clang-format off
  return ((u0=='m' && (!units[1] || strStartsWith(units, "meter")))
       || (!units[2] && units[1]=='m' && (u0=='c' || u0=='d' || u0=='k'))
       || (strStartsWith(units, "decimeter"))
       || (strStartsWith(units, "centimeter"))
       || (strStartsWith(units, "millimeter"))
       || (strStartsWith(units, "kilometer")));
  // clang-format on
}

bool
is_pressure_units(const char *units)
{
  // clang-format off
  return (strStartsWith(units, "millibar")
       || strStartsWith(units, "mb")
       || strStartsWith(units, "hectopas")
       || strStartsWith(units, "hPa")
       || strStartsWith(units, "pa")
       || strStartsWith(units, "Pa"));
  // clang-format on
}

bool
is_DBL_axis(const char *longname)
{
  // clang-format off
  return (str_is_equal(longname, "depth below land")
       || str_is_equal(longname, "depth_below_land")
       || str_is_equal(longname, "levels below the surface"));
  // clang-format on
}

bool
is_depth_axis(const char *stdname, const char *longname)
{
  // clang-format off
  return (str_is_equal(stdname, "depth")
       || str_is_equal(longname, "depth_below_sea")
       || str_is_equal(longname, "depth below sea"));
  // clang-format ofn
}


bool is_height_axis(const char *stdname, const char *longname)
{
  // clang-format off
  return (str_is_equal(stdname, "height")
       || str_is_equal(longname, "height")
       || str_is_equal(longname, "height above the surface"));
  // clang-format on
}

bool
is_altitude_axis(const char *stdname, const char *longname)
{
  // clang-format off
  return (str_is_equal(stdname, "altitude")
       || str_is_equal(longname, "altitude"));
  // clang-format on
}

bool
is_reference_axis(const char *stdname, const char *longname)
{
  // clang-format off
  return ((str_is_equal(longname, "generalized_height") || str_is_equal(longname, "generalized height"))
        && str_is_equal(stdname, "height"));
  // clang-format on
}

bool
is_lon_axis(const char *units, const char *stdname)
{
  bool status = false;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  str_to_lower(lc_units);

  if ((strStartsWith(lc_units, "degree") || strStartsWith(lc_units, "radian"))
      && (strStartsWith(stdname, "grid_longitude") || strStartsWith(stdname, "longitude")))
    {
      status = true;
    }
  else if (strStartsWith(lc_units, "degree") && !strStartsWith(stdname, "grid_latitude") && !strStartsWith(stdname, "latitude"))
    {
      int ioff = 6;
      if (lc_units[ioff] == 's') ioff++;
      if (lc_units[ioff] == ' ') ioff++;
      if (lc_units[ioff] == '_') ioff++;
      if (lc_units[ioff] == 'e') status = true;
    }

  return status;
}

bool
is_lat_axis(const char *units, const char *stdname)
{
  bool status = false;
  char lc_units[16];

  memcpy(lc_units, units, 15);
  lc_units[15] = 0;
  str_to_lower(lc_units);

  if ((strStartsWith(lc_units, "degree") || strStartsWith(lc_units, "radian"))
      && (strStartsWith(stdname, "grid_latitude") || strStartsWith(stdname, "latitude")))
    {
      status = true;
    }
  else if (strStartsWith(lc_units, "degree") && !strStartsWith(stdname, "grid_longitude") && !strStartsWith(stdname, "longitude"))
    {
      int ioff = 6;
      if (lc_units[ioff] == 's') ioff++;
      if (lc_units[ioff] == ' ') ioff++;
      if (lc_units[ioff] == '_') ioff++;
      if (lc_units[ioff] == 'n' || lc_units[ioff] == 's') status = true;
    }

  return status;
}

bool
is_x_axis(const char *units, const char *stdname)
{
  (void) units;
  return (str_is_equal(stdname, "projection_x_coordinate"));
}

bool
is_y_axis(const char *units, const char *stdname)
{
  (void) units;
  return (str_is_equal(stdname, "projection_y_coordinate"));
}

void
cdf_set_gridtype(const char *attstring, int *gridtype)
{
  // clang-format off
  if      (str_is_equal(attstring, "gaussian_reduced")) *gridtype = GRID_GAUSSIAN_REDUCED;
  else if (str_is_equal(attstring, "gaussian"))         *gridtype = GRID_GAUSSIAN;
  else if (strStartsWith(attstring, "spectral"))      *gridtype = GRID_SPECTRAL;
  else if (strStartsWith(attstring, "fourier"))       *gridtype = GRID_FOURIER;
  else if (str_is_equal(attstring, "trajectory"))       *gridtype = GRID_TRAJECTORY;
  else if (str_is_equal(attstring, "generic"))          *gridtype = GRID_GENERIC;
  else if (str_is_equal(attstring, "cell"))             *gridtype = GRID_UNSTRUCTURED;
  else if (str_is_equal(attstring, "unstructured"))     *gridtype = GRID_UNSTRUCTURED;
  else if (str_is_equal(attstring, "curvilinear")) ;
  else if (str_is_equal(attstring, "characterxy"))      *gridtype = GRID_CHARXY;
  else if (str_is_equal(attstring, "sinusoidal")) ;
  else if (str_is_equal(attstring, "laea")) ;
  else if (str_is_equal(attstring, "lcc2")) ;
  else if (str_is_equal(attstring, "linear")) ; // ignore grid type linear
  else
    {
      static bool warn = true;
      if (warn)
        {
          warn = false;
          Warning("NetCDF attribute grid_type='%s' unsupported!", attstring);
        }
    }
  // clang-format on
}

void
cdf_set_zaxistype(const char *attstring, int *zaxistype)
{
  // clang-format off
  if      (str_is_equal(attstring, "toa"))              *zaxistype = ZAXIS_TOA;
  else if (str_is_equal(attstring, "tropopause"))       *zaxistype = ZAXIS_TROPOPAUSE;
  else if (str_is_equal(attstring, "cloudbase"))        *zaxistype = ZAXIS_CLOUD_BASE;
  else if (str_is_equal(attstring, "cloudtop"))         *zaxistype = ZAXIS_CLOUD_TOP;
  else if (str_is_equal(attstring, "isotherm0"))        *zaxistype = ZAXIS_ISOTHERM_ZERO;
  else if (str_is_equal(attstring, "seabottom"))        *zaxistype = ZAXIS_SEA_BOTTOM;
  else if (str_is_equal(attstring, "lakebottom"))       *zaxistype = ZAXIS_LAKE_BOTTOM;
  else if (str_is_equal(attstring, "sedimentbottom"))   *zaxistype = ZAXIS_SEDIMENT_BOTTOM;
  else if (str_is_equal(attstring, "sedimentbottomta")) *zaxistype = ZAXIS_SEDIMENT_BOTTOM_TA;
  else if (str_is_equal(attstring, "sedimentbottomtw")) *zaxistype = ZAXIS_SEDIMENT_BOTTOM_TW;
  else if (str_is_equal(attstring, "mixlayer"))         *zaxistype = ZAXIS_MIX_LAYER;
  else if (str_is_equal(attstring, "atmosphere"))       *zaxistype = ZAXIS_ATMOSPHERE;
  else
    {
      static bool warn = true;
      if (warn)
        {
          warn = false;
          Warning("NetCDF attribute level_type='%s' unsupported!", attstring);
        }
    }
  // clang-format on
}

int
attribute_to_calendar(const char *attstring)
{
  int calendar = CALENDAR_STANDARD;
  // clang-format off
  if      (strStartsWith(attstring, "standard"))  calendar = CALENDAR_STANDARD;
  else if (strStartsWith(attstring, "gregorian")) calendar = CALENDAR_GREGORIAN;
  else if (strStartsWith(attstring, "none"))      calendar = CALENDAR_NONE;
  else if (strStartsWith(attstring, "proleptic")) calendar = CALENDAR_PROLEPTIC;
  else if (strStartsWith(attstring, "360"))       calendar = CALENDAR_360DAYS;
  else if (strStartsWith(attstring, "365") ||
           strStartsWith(attstring, "noleap"))    calendar = CALENDAR_365DAYS;
  else if (strStartsWith(attstring, "366") ||
           strStartsWith(attstring, "all_leap"))  calendar = CALENDAR_366DAYS;
  else Warning("calendar >%s< unsupported!", attstring);
  // clang-format on
  return calendar;
}
