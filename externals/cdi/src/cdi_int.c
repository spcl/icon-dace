#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>

#include "binary.h"
#include "cdf.h"
#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "file.h"
#include "gribapi.h"

#ifdef HAVE_LIBCGRIBEX
#include "cgribex.h"
#endif

int CDI_Default_Calendar = CALENDAR_PROLEPTIC;

int CDI_Default_InstID = CDI_UNDEFID;
int CDI_Default_ModelID = CDI_UNDEFID;
int CDI_Default_TableID = CDI_UNDEFID;
// int cdiNcMissingValue  = CDI_UNDEFID;
int CDI_Netcdf_Chunksizehint = CDI_UNDEFID;
int CDI_Split_Ltype105 = CDI_UNDEFID;

bool CDI_Ignore_Att_Coordinates = false;
bool CDI_Coordinates_Lon_Lat = false;
bool CDI_Ignore_Valid_Range = false;
int CDI_Skip_Records = 0;
const char *CDI_GRIB1_Template = NULL;
const char *CDI_GRIB2_Template = NULL;
int CDI_Convention = CDI_CONVENTION_ECHAM;
int CDI_Inventory_Mode = 1;
int CDI_Version_Info = 1;
int CDI_Query_Abort = 1;
int CDI_Convert_Cubesphere = 1;
int CDI_Read_Cell_Corners = 1;
int CDI_CMOR_Mode = 0;
int CDI_Reduce_Dim = 0;
int CDI_Shuffle = 0;
int CDI_Test = 0;
size_t CDI_Netcdf_Hdr_Pad = 0UL;
size_t CDI_Chunk_Cache = 0UL;
size_t CDI_Chunk_Cache_Max = 0UL;
bool CDI_Netcdf_Lazy_Grid_Load = false;

char *cdiPartabPath = NULL;
int cdiPartabIntern = 1;

double CDI_Default_Missval = -9.E33;
double CDI_Grid_Missval = -9999.;

// clang-format off
static const char Filetypes[][9] = {
  "UNKNOWN",
  "GRIB",
  "GRIB2",
  "NetCDF",
  "NetCDF2",
  "NetCDF4",
  "NetCDF4c",
  "NetCDF5",
  "SERVICE",
  "EXTRA",
  "IEG",
  "NCZarr",
  "HDF5",
};
// clang-format on

int CDI_Debug = 0;  // If set to 1, debugging
int CDI_Recopt = 0;

bool CDI_gribapi_debug = false;
bool CDI_gribapi_grib1 = false;
int cdiDefaultLeveltype = -1;
int cdiDataUnreduced = 0;
int cdiSortName = 0;
int cdiSortParam = 0;
int cdiHaveMissval = 0;

static long
cdi_getenv_int(const char *envName)
{
  long envValue = -1;

  char *envString = getenv(envName);
  if (envString)
    {
      long fact = 1;
      int len = (int) strlen(envString);
      for (int loop = 0; loop < len; loop++)
        {
          if (!isdigit((int) envString[loop]))
            {
              switch (tolower((int) envString[loop]))
                {
                case 'k': fact = 1024; break;
                case 'm': fact = 1048576; break;
                case 'g': fact = 1073741824; break;
                default:
                  fact = 0;
                  Warning("Invalid number string in %s: %s", envName, envString);
                  Warning("%s must comprise only digits [0-9].", envName);
                  break;
                }
              break;
            }
        }

      if (fact) envValue = fact * atol(envString);

      if (CDI_Debug) Message("set %s to %ld", envName, envValue);
    }

  return envValue;
}

static void
cdiPrintDefaults(void)
{
  fprintf(stderr,
          "default instID     :  %d\n"
          "default modelID    :  %d\n"
          "default tableID    :  %d\n"
          "default missval    :  %g\n",
          CDI_Default_InstID, CDI_Default_ModelID, CDI_Default_TableID, CDI_Default_Missval);
}

#ifdef HAVE_LIBFDB5
#include <fdb5/fdb5_config.h>
#endif

void
cdiPrintVersion(void)
{
  fprintf(stdout, "     CDI library version : %s\n", cdiLibraryVersion());
#ifdef HAVE_LIBCGRIBEX
  fprintf(stdout, " cgribex library version : %s\n", cgribexLibraryVersion());
#endif
#ifdef HAVE_LIBGRIB_API
  fprintf(stdout, " ecCodes library version : %s\n", gribapiLibraryVersionString());
#endif
#ifdef HAVE_LIBNETCDF
  fprintf(stdout, "  NetCDF library version : %s\n", cdfLibraryVersion());
#endif
#ifdef HAVE_LIBSERVICE
  fprintf(stdout, "    exse library version : %s\n", srvLibraryVersion());
#endif
  fprintf(stdout, "    FILE library version : %s\n", fileLibraryVersion());
#ifdef HAVE_LIBFDB5
  fprintf(stdout, "    FDB5 library version : %s\n", fdb5_version());
#endif
}

static void
cdiPrintDatatypes(void)
{
#define XSTRING(x) #x
#define STRING(x) XSTRING(x)

  fprintf(stderr,
          "+-------------+-------+\n"
          "| types       | bytes |\n"
          "+-------------+-------+\n"
          "| void *      |   %3d |\n"
          "+-------------+-------+\n"
          "| char        |   %3d |\n"
          "+-------------+-------+\n"
          "| bool        |   %3d |\n"
          "| short       |   %3d |\n"
          "| int         |   %3d |\n"
          "| long        |   %3d |\n"
          "| long long   |   %3d |\n"
          "| size_t      |   %3d |\n"
          "| off_t       |   %3d |\n"
          "+-------------+-------+\n"
          "| float       |   %3d |\n"
          "| double      |   %3d |\n"
          "| long double |   %3d |\n"
          "+-------------+-------+\n\n",
          (int) sizeof(void *), (int) sizeof(char), (int) sizeof(bool), (int) sizeof(short), (int) sizeof(int), (int) sizeof(long),
          (int) sizeof(long long), (int) sizeof(size_t), (int) sizeof(off_t), (int) sizeof(float), (int) sizeof(double),
          (int) sizeof(long double));

  fprintf(stderr,
          "+-------------+-----------+\n"
          "| INT32       | %-9s |\n"
          "| INT64       | %-9s |\n"
          "| FLT32       | %-9s |\n"
          "| FLT64       | %-9s |\n"
          "| SizeType    | %-9s |\n"
          "+-------------+-----------+\n",
          STRING(INT32), STRING(INT64), STRING(FLT32), STRING(FLT64), STRING(CDI_SIZE_TYPE));

  fprintf(stderr, "\n  byte ordering is %s\n\n",
          ((HOST_ENDIANNESS == CDI_BIGENDIAN)
               ? "BIGENDIAN"
               : ((HOST_ENDIANNESS == CDI_LITTLEENDIAN) ? "LITTLEENDIAN" : "Unhandled endianness!")));

#undef STRING
#undef XSTRING
}

void
cdiDebug(int level)
{
  unsigned ulevel = (level == 1) ? (1U << 16) : (unsigned) level;

  if (ulevel & 2) CDI_Debug = 1;
  if (ulevel & 4) memDebug(1);
  if (ulevel & 8) fileDebug(1);

  if (ulevel & 16)
    {
#ifdef HAVE_LIBCGRIBEX
      gribSetDebug(1);
#endif
#ifdef HAVE_LIBNETCDF
      cdfDebug(1);
#endif
#ifdef HAVE_LIBSERVICE
      srvDebug(1);
#endif
#ifdef HAVE_LIBEXTRA
      extDebug(1);
#endif
#ifdef HAVE_LIBIEG
      iegDebug(1);
#endif
    }

  if (CDI_Debug)
    {
      cdiPrintDefaults();
      cdiPrintDatatypes();
    }
}

int
cdiHaveFiletype(int filetype)
{
  int status = 0;

  switch (filetype)
    {
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: status = 1; break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: status = 1; break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: status = 1; break;
#endif
#ifdef HAVE_LIBGRIB
#if defined HAVE_LIBGRIB_API || defined HAVE_LIBCGRIBEX
    case CDI_FILETYPE_GRB: status = 1; break;
#endif
#ifdef HAVE_LIBGRIB_API
    case CDI_FILETYPE_GRB2: status = 1; break;
#endif
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC: status = 1; break;
#ifdef HAVE_NETCDF2
    case CDI_FILETYPE_NC2: status = 1; break;
#endif
#ifdef HAVE_NETCDF4
    case CDI_FILETYPE_NC4: status = 1; break;
    case CDI_FILETYPE_NC4C: status = 1; break;
#endif
#ifdef HAVE_NCZARR
    case CDI_FILETYPE_NCZARR: status = 1; break;
#endif
#ifdef HAVE_NETCDF5
    case CDI_FILETYPE_NC5: status = 1; break;
#endif
#endif
    default: status = 0; break;
    }

  return status;
}

void
cdiDefTableID(int tableID)
{
  CDI_Default_TableID = tableID;
  int modelID = CDI_Default_ModelID = tableInqModel(tableID);
  CDI_Default_InstID = modelInqInstitut(modelID);
}

void
cdiSetEccodesGrib1(bool value)
{
#ifndef HAVE_LIBGRIB_API
  if (value)
    {
      Warning("ecCodes support not compiled in, used CGRIBEX to decode/encode GRIB1 records!");
      value = false;
    }
#endif
  CDI_gribapi_grib1 = value;
}

void
cdiInitialize(void)
{
  static bool Init_CDI = false;

  if (!Init_CDI)
    {
      Init_CDI = true;
      char *envstr;
      long value;

#ifdef HAVE_LIBCGRIBEX
      gribFixZSE(1);    // 1: Fix ZeroShiftError of simple packed spherical harmonics
      gribSetConst(1);  // 1: Don't pack constant fields on regular grids
#endif
#ifdef HAVE_LIBGRIB_API
      grib_multi_support_off(NULL);
#endif

      value = cdi_getenv_int("CDI_DEBUG");
      if (value >= 0) CDI_Debug = (int) value;

      value = cdi_getenv_int("CDI_GRIBAPI_DEBUG");
      if (value >= 0) CDI_gribapi_debug = (bool) value;

      value = cdi_getenv_int("CDI_ECCODES_DEBUG");
      if (value >= 0) CDI_gribapi_debug = (bool) value;

      value = cdi_getenv_int("CDI_ECCODES_GRIB1");
      if (value >= 0) cdiSetEccodesGrib1((bool) value);

      value = cdi_getenv_int("CDI_READ_CELL_CORNERS");
      if (value >= 0) CDI_Read_Cell_Corners = (int) value;

      value = cdi_getenv_int("CDI_RECOPT");
      if (value >= 0) CDI_Recopt = (int) value;

      value = cdi_getenv_int("CDI_REGULARGRID");
      if (value >= 0) cdiDataUnreduced = (int) value;

      value = cdi_getenv_int("CDI_SORTNAME");
      if (value >= 0) cdiSortName = (int) value;

      value = cdi_getenv_int("CDI_SORTPARAM");
      if (value >= 0) cdiSortParam = (int) value;

      value = cdi_getenv_int("CDI_HAVE_MISSVAL");
      if (value >= 0) cdiHaveMissval = (int) value;

      value = cdi_getenv_int("CDI_LEVELTYPE");
      if (value >= 0) cdiDefaultLeveltype = (int) value;

      value = cdi_getenv_int("CDI_NETCDF_HDR_PAD");
      if (value >= 0) CDI_Netcdf_Hdr_Pad = (size_t) value;

      value = cdi_getenv_int("CDI_CHUNK_CACHE");
      if (value >= 0) CDI_Chunk_Cache = (size_t) value;

      value = cdi_getenv_int("CDI_CHUNK_CACHE_MAX");
      if (value >= 0) CDI_Chunk_Cache_Max = (size_t) value;

      value = cdi_getenv_int("CDI_TEST");
      if (value >= 0) CDI_Test = (size_t) value;

      envstr = getenv("CDI_GRIB1_TEMPLATE");
      if (envstr) CDI_GRIB1_Template = envstr;

      envstr = getenv("CDI_GRIB2_TEMPLATE");
      if (envstr) CDI_GRIB2_Template = envstr;

      envstr = getenv("CDI_SHUFFLE");
      if (envstr) CDI_Shuffle = atoi(envstr);

      envstr = getenv("CDI_MISSVAL");
      if (envstr) CDI_Default_Missval = atof(envstr);
      /*
      envstr = getenv("NC_MISSING_VALUE");
      if ( envstr ) cdiNcMissingValue = atoi(envstr);
      */
      envstr = getenv("NC_CHUNKSIZEHINT");
      if (envstr) CDI_Netcdf_Chunksizehint = atoi(envstr);

      envstr = getenv("SPLIT_LTYPE_105");
      if (envstr) CDI_Split_Ltype105 = atoi(envstr);

      envstr = getenv("IGNORE_ATT_COORDINATES");
      if (envstr) CDI_Ignore_Att_Coordinates = atoi(envstr) > 0;

      envstr = getenv("CDI_COORDINATES_LONLAT");
      if (envstr) CDI_Coordinates_Lon_Lat = atoi(envstr) > 0;

      envstr = getenv("IGNORE_VALID_RANGE");
      if (envstr) CDI_Ignore_Valid_Range = atoi(envstr) > 0;

      envstr = getenv("CDI_SKIP_RECORDS");
      if (envstr)
        {
          CDI_Skip_Records = atoi(envstr);
          CDI_Skip_Records = CDI_Skip_Records > 0 ? CDI_Skip_Records : 0;
        }

      envstr = getenv("CDI_CONVENTION");
      if (envstr)
        {
          if (str_is_equal(envstr, "CF") || str_is_equal(envstr, "cf"))
            {
              CDI_Convention = CDI_CONVENTION_CF;
              if (CDI_Debug) Message("CDI convention was set to CF!");
            }
        }

      envstr = getenv("CDI_INVENTORY_MODE");
      if (envstr)
        {
          if (strncmp(envstr, "time", 4) == 0)
            {
              CDI_Inventory_Mode = 2;
              if (CDI_Debug) Message("Inventory mode was set to timestep!");
            }
        }

      envstr = getenv("CDI_QUERY_ABORT");
      if (envstr)
        {
          const int ival = atoi(envstr);
          if (ival == 0 || ival == 1)
            {
              CDI_Query_Abort = ival;
              if (CDI_Debug) Message("CDI_Query_Abort = %s", envstr);
            }
        }

      envstr = getenv("CDI_VERSION_INFO");
      if (envstr)
        {
          const int ival = atoi(envstr);
          if (ival == 0 || ival == 1)
            {
              CDI_Version_Info = ival;
              if (CDI_Debug) Message("CDI_Version_Info = %s", envstr);
            }
        }

      envstr = getenv("CDI_CONVERT_CUBESPHERE");
      if (envstr)
        {
          const int ival = atoi(envstr);
          if (ival == 0 || ival == 1)
            {
              CDI_Convert_Cubesphere = ival;
              if (CDI_Debug) Message("CDI_Convert_Cubesphere = %s", envstr);
            }
        }

      envstr = getenv("CDI_CALENDAR");
      if (envstr)
        {
          // clang-format off
	  if      (strncmp(envstr, "standard", 8)  == 0) CDI_Default_Calendar = CALENDAR_STANDARD;
	  else if (strncmp(envstr, "gregorian", 9) == 0) CDI_Default_Calendar = CALENDAR_GREGORIAN;
	  else if (strncmp(envstr, "proleptic", 9) == 0) CDI_Default_Calendar = CALENDAR_PROLEPTIC;
	  else if (strncmp(envstr, "360days", 7)   == 0) CDI_Default_Calendar = CALENDAR_360DAYS;
	  else if (strncmp(envstr, "365days", 7)   == 0) CDI_Default_Calendar = CALENDAR_365DAYS;
	  else if (strncmp(envstr, "366days", 7)   == 0) CDI_Default_Calendar = CALENDAR_366DAYS;
	  else if (strncmp(envstr, "none", 4)      == 0) CDI_Default_Calendar = CALENDAR_NONE;
          // clang-format on
          if (CDI_Debug) Message("Default calendar set to %s!", envstr);
        }
#ifdef HAVE_LIBCGRIBEX
      gribSetCalendar(CDI_Default_Calendar);
#endif

      envstr = getenv("PARTAB_INTERN");
      if (envstr) cdiPartabIntern = atoi(envstr);

      envstr = getenv("PARTAB_PATH");
      if (envstr) cdiPartabPath = strdup(envstr);
    }
}

const char *
strfiletype(int filetype)
{
  const int size = (int) (sizeof(Filetypes) / sizeof(char *));
  return (filetype > 0 && filetype < size) ? Filetypes[filetype] : Filetypes[0];
}

void
cdiDefGlobal(const char *string, int value)
{
  // clang-format off
  if      (str_is_equal(string, "REGULARGRID")          ) cdiDataUnreduced = value;
  else if (str_is_equal(string, "ECCODES_DEBUG")        ) CDI_gribapi_debug = (bool) value;
  else if (str_is_equal(string, "ECCODES_GRIB1")        ) cdiSetEccodesGrib1((bool) value);
  else if (str_is_equal(string, "SORTNAME")             ) cdiSortName = value;
  else if (str_is_equal(string, "SORTPARAM")            ) cdiSortParam = value;
  else if (str_is_equal(string, "HAVE_MISSVAL")         ) cdiHaveMissval = value;
  else if (str_is_equal(string, "NC_CHUNKSIZEHINT")     ) CDI_Netcdf_Chunksizehint = value;
  else if (str_is_equal(string, "READ_CELL_CORNERS")    ) CDI_Read_Cell_Corners = value;
  else if (str_is_equal(string, "CMOR_MODE")            ) CDI_CMOR_Mode = value;
  else if (str_is_equal(string, "REDUCE_DIM")           ) CDI_Reduce_Dim = value;
  else if (str_is_equal(string, "NETCDF_HDR_PAD")       ) CDI_Netcdf_Hdr_Pad = (size_t) value;
  else if (str_is_equal(string, "NETCDF_LAZY_GRID_LOAD")) CDI_Netcdf_Lazy_Grid_Load = (bool) value;
  else Warning("Unsupported global key: %s", string);
  // clang-format on
}

void
cdiDefMissval(double missval)
{
  cdiInitialize();

  CDI_Default_Missval = missval;
}

double
cdiInqMissval(void)
{
  cdiInitialize();

  return CDI_Default_Missval;
}

bool
cdiFiletypeIsExse(int filetype)
{
  return (filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG);
}

int
cdiBaseFiletype(int filetype)
{
  switch (filetype)
    {
    case CDI_FILETYPE_GRB:
    case CDI_FILETYPE_GRB2: return CDI_FILETYPE_GRIB;
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NCZARR:
    case CDI_FILETYPE_NC5: return CDI_FILETYPE_NETCDF;
    }

  return filetype;
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
