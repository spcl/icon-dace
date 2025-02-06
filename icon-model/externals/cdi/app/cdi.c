#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include <cdi.h>

int vlistInqVarMissvalUsed(int vlistID, int varID);
#ifndef DBL_IS_NAN
#if defined(HAVE_DECL_ISNAN)
#define DBL_IS_NAN(x) (isnan(x))
#elif defined(FP_NAN)
#define DBL_IS_NAN(x) (fpclassify(x) == FP_NAN)
#else
#define DBL_IS_NAN(x) ((x) != (x))
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#define DBL_IS_EQUAL(x, y) (DBL_IS_NAN(x) || DBL_IS_NAN(y) ? (DBL_IS_NAN(x) && DBL_IS_NAN(y) ? 1 : 0) : !(x < y || y < x))
#endif

#ifndef IS_NOT_EQUAL
#define IS_NOT_EQUAL(x, y) (x < y || y < x)
#endif

#ifndef HOST_ENDIANNESS
#ifdef __cplusplus
static const uint32_t HOST_ENDIANNESS_temp[1] = { UINT32_C(0x00030201) };
#define HOST_ENDIANNESS (((const unsigned char *) HOST_ENDIANNESS_temp)[0])
#else
#define HOST_ENDIANNESS (((const unsigned char *) &(const uint32_t[1]){ UINT32_C(0x00030201) })[0])
#endif
#endif

#include "printinfo.h"

#ifdef __cplusplus
extern "C"
{
#endif
  void cdiDefTableID(int tableID);
#ifdef __cplusplus
}
#endif

int getopt(int argc, char *const argv[], const char *optstring);

extern char *optarg;
extern int optind, opterr, optopt;

static char *Progname;

static int DefaultFileType = CDI_UNDEFID;
static int DefaultDataType = CDI_UNDEFID;
static int DefaultByteorder = CDI_UNDEFID;

static int filterId = 0;
static int params[8];
static int maxParams = sizeof(params) / sizeof(params[0]);
static int numParams = 0;
static int comptype = CDI_COMPRESS_NONE;  // Compression type
static int complevel = 0;                 // Compression level

enum datamode
{
  SP_MODE,
  DP_MODE
};
static int datamode = DP_MODE;

static void
version(void)
{
  int filetypes[] = { CDI_FILETYPE_SRV, CDI_FILETYPE_EXT, CDI_FILETYPE_IEG,  CDI_FILETYPE_GRB, CDI_FILETYPE_GRB2,  CDI_FILETYPE_NC,
                      CDI_FILETYPE_NC2, CDI_FILETYPE_NC4, CDI_FILETYPE_NC4C, CDI_FILETYPE_NC5, CDI_FILETYPE_NCZARR };
  const char *typenames[] = { "srv", "ext", "ieg", "grb", "grb2", "nc", "nc2", "nc4", "nc4c", "nc5", "nczarr" };

  fprintf(stderr, "CDI version 2.2\n");
#ifdef COMPILER
  fprintf(stderr, "C Compiler: %s\n", COMPILER);
#endif
#ifdef COMP_VERSION
  fprintf(stderr, "C version: %s\n", COMP_VERSION);
#endif
#ifdef SYSTEM_TYPE
  fprintf(stderr, "System: %s\n", SYSTEM_TYPE);
#endif

  fprintf(stderr, "Features:");
#ifdef HAVE_CF_INTERFACE
  fprintf(stderr, " Fortran");
#endif
#ifdef HAVE_LIBPTHREAD
  fprintf(stderr, " PTHREADS");
#endif
#ifdef _OPENMP
  fprintf(stderr, " OpenMP");
#endif
#ifdef HAVE_NETCDF4
  fprintf(stderr, " NC4");
#ifdef HAVE_NC4HDF5
  fprintf(stderr, "/HDF5");
#ifdef HAVE_NC4HDF5_THREADSAFE
  fprintf(stderr, "/threadsafe");
#endif
#endif
#endif
#ifdef HAVE_LIBNC_DAP
  fprintf(stderr, " OPeNDAP");
#endif
#ifdef HAVE_LIBSZ
  fprintf(stderr, " SZ");
#endif
#ifdef HAVE_LIBJASPER
  fprintf(stderr, " JASPER");
#endif
#ifdef HAVE_LIBDRMAA
  fprintf(stderr, " DRMAA");
#endif
#ifdef HAVE_LIBCURL
  fprintf(stderr, " CURL");
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "Filetypes: ");
  for (size_t i = 0; i < sizeof(filetypes) / sizeof(int); ++i)
    if (cdiHaveFiletype(filetypes[i])) fprintf(stderr, "%s ", typenames[i]);
  fprintf(stderr, "\n");

  cdiPrintVersion();
  fprintf(stderr, "\n");
  /*
    1.0.0   6 Feb 2001 : initial version
    1.1.0  30 Jul 2003 : missing values implemented
    1.2.0   8 Aug 2003 : changes for CDI library version 0.7.0
    1.3.0  10 Feb 2004 : changes for CDI library version 0.7.9
    1.4.0   5 May 2004 : changes for CDI library version 0.8.1 (error handling)
    1.4.1  18 Sep 2004 : netCDF 2 support
    1.4.2  22 Mar 2005 : change level from int to double
    1.4.3  11 Apr 2005 : change date and time format to ISO
    1.5.0  22 Nov 2005 : IEG support
    1.5.1  21 Feb 2006 : add option -s for short info
    1.6.0   1 Aug 2006 : add option -z szip for SZIP compression of GRIB records
    1.6.1  27 Feb 2007 : short info with ltype for GENERIC zaxis
    1.6.2   3 Jan 2008 : changes for CDI library version 1.1.0 (compress)
    1.6.3  26 Mar 2008 : call streamDefTimestep also if ntsteps = 0 (buf fix)
    1.7.0  11 Apr 2008 : add option -z zip for deflate compression of netCDF4 variables
    1.7.1   1 Nov 2009 : add option -z jpeg for JPEG compression of GRIB2 records
    1.7.2  14 Nov 2012 : add optional compression level -z zip[_1-9]
    1.9.0  29 May 2019 : add option -i to set number of input worker
    2.0.0  02 Feb 2022 : changed date/time handling to CdiDateTime
    2.1.0  17 Jun 2022 : add NCZarr support
    2.2.0  23 Jul 2023 : add option -E and -T
  */
}

static void
usage(void)
{
  const char *name;

  fprintf(stderr, "usage : %s  [Option]  [ifile]  [ofile]\n", Progname);

  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  fprintf(stderr, "    -d             Print debugging information\n");
  fprintf(stderr, "    -E             Use ecCodes to decode/encode GRIB1 messages\n");
  fprintf(stderr, "    -f <format>    Format of the output file. (grb, grb2, nc, nc2, nc4, nc4c, nc5, nczarr, srv, ext or ieg)\n");
  fprintf(stderr, "    -i <num>       Number of worker to decode/decompress GRIB records\n");
  fprintf(stderr, "    -m             Move records\n");
  fprintf(stderr, "    -r             Use CDI record API\n");
  fprintf(stderr, "    -s             give short information if ofile is missing\n");
  fprintf(stderr, "    -T             Pre scan hole GRIB file to get the number of timesteps\n");
  fprintf(stderr, "    -t <table>     Parameter table name/file\n");
  fprintf(stderr, "                   Predefined tables: ");
  for (int id = 0; id < tableInqNumber(); id++)
    if ((name = tableInqNamePtr(id))) fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    -V             Print version number\n");
  fprintf(stderr, "    -z szip        SZIP compression of GRIB1 records\n");
  fprintf(stderr, "       jpeg        JPEG compression of GRIB2 records\n");
  fprintf(stderr, "        zip[_1-9]  Deflate compression of netCDF4 variables\n");
  fprintf(stderr, "        zstd[_1-19] zstandard compression of netCDF4 variables\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Report bugs to <https://code.mpimet.mpg.de/projects/cdi>\n");
}

static void
printInfo(CdiDateTime vdatetime, char *varname, double level, size_t datasize, int number, size_t nmiss, double missval,
          const double *data, int vardis)
{
  static int rec = 0;
  size_t ivals = 0, imiss = 0;
  double arrmean, arrmin, arrmax;

  if (!rec)
    {
      if (vardis)
        fprintf(stdout,
                "   Rec :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name\n");
      //         ----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+
      else
        fprintf(stdout,
                "   Rec :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter ID\n");
      //         ----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+
    }

  char vdatestr[32], vtimestr[32];
  date2str(vdatetime.date, vdatestr, sizeof(vdatestr));
  time2str(vdatetime.time, vtimestr, sizeof(vtimestr));

  fprintf(stdout, "%6d :%s %s %7g ", ++rec, vdatestr, vtimestr, level);
  fprintf(stdout, "%8zu ", datasize);
  fprintf(stdout, "%7zu :", nmiss);

  if (number == CDI_REAL)
    {
      if (nmiss > 0)
        {
          arrmean = 0;
          arrmin = 1.e300;
          arrmax = -1.e300;
          for (size_t i = 0; i < datasize; i++)
            {
              if (!DBL_IS_EQUAL(data[i], missval))
                {
                  if (data[i] < arrmin) arrmin = data[i];
                  if (data[i] > arrmax) arrmax = data[i];
                  arrmean += data[i];
                  ivals++;
                }
            }
          imiss = datasize - ivals;
          datasize = ivals;
        }
      else
        {
          arrmean = data[0];
          arrmin = data[0];
          arrmax = data[0];
          for (size_t i = 1; i < datasize; i++)
            {
              if (data[i] < arrmin) arrmin = data[i];
              if (data[i] > arrmax) arrmax = data[i];
              arrmean += data[i];
            }
        }

      if (datasize > 0) arrmean /= datasize;

      fprintf(stdout, "%#12.5g%#12.5g%#12.5g", arrmin, arrmean, arrmax);
    }
  else
    {
      size_t nvals_r = 0, nvals_i = 0;
      double arrsum_r = 0, arrsum_i = 0, arrmean_r = 0, arrmean_i = 0;

      for (size_t i = 0; i < datasize; i++)
        {
          if (!DBL_IS_EQUAL(data[i * 2], missval))
            {
              arrsum_r += data[i * 2];
              nvals_r++;
            }
          if (!DBL_IS_EQUAL(data[i * 2 + 1], missval))
            {
              arrsum_i += data[i * 2 + 1];
              nvals_i++;
            }
        }

      imiss = datasize - nvals_r;

      if (nvals_r > 0) arrmean_r = arrsum_r / nvals_r;
      if (nvals_i > 0) arrmean_i = arrsum_i / nvals_i;
      fprintf(stdout, "  -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
    }

  fprintf(stdout, " : %-14s\n", varname);

  if (imiss != nmiss && nmiss > 0) fprintf(stdout, "Found %zu of %zu missing values!\n", imiss, nmiss);
}

static const char *
tunit2str(int tunits)
{
  // clang-format off
  if      (tunits == TUNIT_YEAR)       return ("years");
  else if (tunits == TUNIT_MONTH)      return ("months");
  else if (tunits == TUNIT_DAY)        return ("days");
  else if (tunits == TUNIT_12HOURS)    return ("12hours");
  else if (tunits == TUNIT_6HOURS)     return ("6hours");
  else if (tunits == TUNIT_3HOURS)     return ("3hours");
  else if (tunits == TUNIT_HOUR)       return ("hours");
  else if (tunits == TUNIT_30MINUTES)  return ("30minutes");
  else if (tunits == TUNIT_QUARTER)    return ("15minutes");
  else if (tunits == TUNIT_MINUTE)     return ("minutes");
  else if (tunits == TUNIT_SECOND)     return ("seconds");
  else                                 return ("unknown");
  // clang-format on
}

static const char *
calendar2str(int calendar)
{
  // clang-format off
  if      (calendar == CALENDAR_STANDARD)  return ("standard");
  else if (calendar == CALENDAR_PROLEPTIC) return ("proleptic_gregorian");
  else if (calendar == CALENDAR_360DAYS)   return ("360_day");
  else if (calendar == CALENDAR_365DAYS)   return ("365_day");
  else if (calendar == CALENDAR_366DAYS)   return ("366_day");
  else                                     return ("unknown");
  // clang-format on
}

static void
limit_string_length(char *string, size_t maxlen)
{
  string[maxlen - 1] = 0;
  size_t len = strlen(string);

  if (len > 10)
    {
      for (size_t i = 3; i < len; ++i)
        if (string[i] == ' ' || string[i] == ',' || (i > 10 && string[i] == '.'))
          {
            string[i] = 0;
            break;
          }
    }
}

static void
print_short_info(int streamID, int vlistID, int vardis)
{
  char tmpname[CDI_MAX_NAME];
  char varname[CDI_MAX_NAME];
  char pstr[4];
  char paramstr[32];

  fprintf(stdout, "   File format");
  fprintf(stdout, " : ");
  printFiletype(streamID, vlistID);

  // vlistPrint(vlistID);
  int nvars = vlistNvars(vlistID);
  int nsubtypes = vlistNsubtypes(vlistID);

  if (nsubtypes > 0)
    fprintf(stdout, "   Var : Institut Source   T Steptype Subtypes Levels Num    Points Num Dtype : ");
  else
    fprintf(stdout, "   Var : Institut Source   T Steptype Levels Num    Points Num Dtype : ");

  if (vardis)
    fprintf(stdout, "Parameter name\n");
  else
    fprintf(stdout, "Parameter ID\n");

  for (int varID = 0; varID < nvars; varID++)
    {
      int param = vlistInqVarParam(vlistID, varID);
      int gridID = vlistInqVarGrid(vlistID, varID);
      int zaxisID = vlistInqVarZaxis(vlistID, varID);

      fprintf(stdout, "%6d : ", varID + 1);

      // institute info
      const char *instptr = institutInqNamePtr(vlistInqVarInstitut(vlistID, varID));
      strcpy(tmpname, "unknown");
      if (instptr) strncpy(tmpname, instptr, CDI_MAX_NAME - 1);
      limit_string_length(tmpname, CDI_MAX_NAME);
      fprintf(stdout, "%-8s ", tmpname);

      // source info
      const char *modelptr = modelInqNamePtr(vlistInqVarModel(vlistID, varID));
      strcpy(tmpname, "unknown");
      if (modelptr) strncpy(tmpname, modelptr, CDI_MAX_NAME - 1);
      limit_string_length(tmpname, CDI_MAX_NAME);
      fprintf(stdout, "%-8s ", tmpname);

      // timetype
      int timetype = vlistInqVarTimetype(vlistID, varID);
      fprintf(stdout, "%c ", timetype == TIME_CONSTANT ? 'c' : 'v');

      // tsteptype
      int tsteptype = vlistInqVarTsteptype(vlistID, varID);
      // clang-format off
	  if      (tsteptype == TSTEP_INSTANT ) fprintf(stdout, "%-8s ", "instant");
	  else if (tsteptype == TSTEP_INSTANT2) fprintf(stdout, "%-8s ", "instant");
	  else if (tsteptype == TSTEP_INSTANT3) fprintf(stdout, "%-8s ", "instant");
	  else if (tsteptype == TSTEP_MIN     ) fprintf(stdout, "%-8s ", "min");
	  else if (tsteptype == TSTEP_MAX     ) fprintf(stdout, "%-8s ", "max");
	  else if (tsteptype == TSTEP_AVG     ) fprintf(stdout, "%-8s ", "avg");
	  else if (tsteptype == TSTEP_ACCUM   ) fprintf(stdout, "%-8s ", "accum");
	  else if (tsteptype == TSTEP_RANGE   ) fprintf(stdout, "%-8s ", "range");
	  else if (tsteptype == TSTEP_DIFF    ) fprintf(stdout, "%-8s ", "diff");
	  else if (tsteptype == TSTEP_SUM     ) fprintf(stdout, "%-8s ", "sum");
	  else                                  fprintf(stdout, "%-8s ", "unknown");
      // clang-format on

      if (nsubtypes > 0)
        {
          int subtypeID = vlistInqVarSubtype(vlistID, varID);
          int subtypesize = subtypeInqSize(subtypeID);
          fprintf(stdout, " %6d  ", subtypesize);
          fprintf(stdout, "%3d ", vlistSubtypeIndex(vlistID, subtypeID) + 1);
        }

      // layer info
      int levelsize = zaxisInqSize(zaxisID);
      fprintf(stdout, "%6d ", levelsize);
      fprintf(stdout, "%3d ", vlistZaxisIndex(vlistID, zaxisID) + 1);

      // grid info
      size_t gridsize = gridInqSize(gridID);
      fprintf(stdout, "%9zu ", gridsize);
      fprintf(stdout, "%3d ", vlistGridIndex(vlistID, gridID) + 1);

      // datatype
      int datatype = vlistInqVarDatatype(vlistID, varID);
      // clang-format off
      if      (datatype == CDI_DATATYPE_PACK  ) strcpy(pstr, "P0");
      else if (datatype > 0 && datatype <= 32 ) snprintf(pstr, sizeof(pstr), "P%d", datatype);
      else if (datatype == CDI_DATATYPE_CPX32 ) strcpy(pstr, "C32");
      else if (datatype == CDI_DATATYPE_CPX64 ) strcpy(pstr, "C64");
      else if (datatype == CDI_DATATYPE_FLT32 ) strcpy(pstr, "F32");
      else if (datatype == CDI_DATATYPE_FLT64 ) strcpy(pstr, "F64");
      else if (datatype == CDI_DATATYPE_INT8  ) strcpy(pstr, "I8");
      else if (datatype == CDI_DATATYPE_INT16 ) strcpy(pstr, "I16");
      else if (datatype == CDI_DATATYPE_INT32 ) strcpy(pstr, "I32");
      else if (datatype == CDI_DATATYPE_UINT8 ) strcpy(pstr, "U8");
      else if (datatype == CDI_DATATYPE_UINT16) strcpy(pstr, "U16");
      else if (datatype == CDI_DATATYPE_UINT32) strcpy(pstr, "U32");
      else                                      strcpy(pstr, "-1");
      // clang-format on

      fprintf(stdout, " %-3s", pstr);

      int compType = vlistInqVarCompType(vlistID, varID);
      bool isCompressed = (compType != CDI_COMPRESS_NONE);
      fprintf(stdout, "%c ", isCompressed ? (int) comptype_to_name(compType)[0] : ' ');

      // parameter info
      fprintf(stdout, ": ");

      cdiParamToString(param, paramstr, sizeof(paramstr));

      if (vardis)
        {
          vlistInqVarName(vlistID, varID, varname);
          fprintf(stdout, "%-14s", varname);
        }
      else
        fprintf(stdout, "%-14s", paramstr);

      fprintf(stdout, "\n");
    }

  fprintf(stdout, "   Grid coordinates");
  fprintf(stdout, " :\n");

  printGridInfo(vlistID);

  fprintf(stdout, "   Vertical coordinates");
  fprintf(stdout, " :\n");

  printZaxisInfo(vlistID);

  if (nsubtypes > 0)
    {
      fprintf(stdout, "   Subtypes");
      fprintf(stdout, " :\n");

      printSubtypeInfo(vlistID);
    }

  int taxisID = vlistInqTaxis(vlistID);
  int ntsteps = vlistNtsteps(vlistID);

  if (ntsteps != 0)
    {
      if (ntsteps == CDI_UNDEFID)
        fprintf(stdout, "   Time coordinate :  unlimited steps\n");
      else
        fprintf(stdout, "   Time coordinate :  %d step%s\n", ntsteps, ntsteps == 1 ? "" : "s");

      if (taxisID != CDI_UNDEFID)
        {
          if (taxisInqType(taxisID) == TAXIS_RELATIVE)
            {
              CdiDateTime dt = taxisInqRdatetime(taxisID);

              fprintf(stdout, "     RefTime = %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d", dt.date.year, dt.date.month, dt.date.day,
                      dt.time.hour, dt.time.minute, dt.time.second);
              if (dt.time.ms) fprintf(stdout, ".%d", dt.time.ms);

              int tunits = taxisInqTunit(taxisID);
              if (tunits != CDI_UNDEFID) fprintf(stdout, "  Units = %s", tunit2str(tunits));

              int calendar = taxisInqCalendar(taxisID);
              if (calendar != CDI_UNDEFID) fprintf(stdout, "  Calendar = %s", calendar2str(calendar));

              if (taxisHasBounds(taxisID)) fprintf(stdout, "  Bounds = true");

              fprintf(stdout, "\n");
            }
        }

      fprintf(stdout, "  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss  YYYY-MM-DD hh:mm:ss\n");

      printTimesteps(streamID, taxisID, 0);

      fprintf(stdout, "\n");
    }
}

static void
setDefaultDataType(char *datatypestr)
{
  int nbits = -1;
  enum
  {
    D_UINT,
    D_INT,
    D_FLT,
    D_CPX
  };
  int dtype = -1;

  int datatype = tolower(*datatypestr);
  // clang-format off
  if      (datatype == 'i') { dtype = D_INT;  datatypestr++; }
  else if (datatype == 'u') { dtype = D_UINT; datatypestr++; }
  else if (datatype == 'f') { dtype = D_FLT;  datatypestr++; }
  else if (datatype == 'c') { dtype = D_CPX;  datatypestr++; }
  else if (datatype == 'p') {                 datatypestr++; }
  // clang-format on

  if (isdigit((int) *datatypestr))
    {
      nbits = atoi(datatypestr);
      datatypestr += 1;
      if (nbits >= 10) datatypestr += 1;

      if (dtype == -1)
        {
          if (nbits > 0 && nbits < 32)
            DefaultDataType = nbits;
          else if (nbits == 32)
            DefaultDataType = (DefaultFileType == CDI_FILETYPE_GRB) ? CDI_DATATYPE_PACK32 : CDI_DATATYPE_FLT32;
          else if (nbits == 64)
            DefaultDataType = CDI_DATATYPE_FLT64;
          else
            {
              fprintf(stderr, "Unsupported number of bits %d!\n", nbits);
              fprintf(
                  stderr,
                  "Use I8/I16/I32/F32/F64 for nc/nc2/nc4/nc4c/nc5/nczarr; F32/F64 for grb2/srv/ext/ieg; P1 - P24 for grb/grb2.\n");
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (dtype == D_INT)
            {
              if (nbits == 8)
                DefaultDataType = CDI_DATATYPE_INT8;
              else if (nbits == 16)
                DefaultDataType = CDI_DATATYPE_INT16;
              else if (nbits == 32)
                DefaultDataType = CDI_DATATYPE_INT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype INT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          /*
          else if ( dtype == D_UINT )
            {
              if      ( nbits ==  8 ) DefaultDataType = CDI_DATATYPE_UINT8;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype UINT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          */
          else if (dtype == D_FLT)
            {
              if (nbits == 32)
                DefaultDataType = CDI_DATATYPE_FLT32;
              else if (nbits == 64)
                DefaultDataType = CDI_DATATYPE_FLT64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype FLT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if (dtype == D_CPX)
            {
              if (nbits == 32)
                DefaultDataType = CDI_DATATYPE_CPX32;
              else if (nbits == 64)
                DefaultDataType = CDI_DATATYPE_CPX64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype CPX!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
        }
    }

  if (*datatypestr != 0)
    {
      if (*datatypestr == 'l' || *datatypestr == 'L')
        {
          if (HOST_ENDIANNESS == CDI_BIGENDIAN) DefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;
        }
      else if (*datatypestr == 'b' || *datatypestr == 'B')
        {
          if (HOST_ENDIANNESS == CDI_LITTLEENDIAN) DefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: >%s< !\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}

static void
setDefaultFileType(char *filetypestr)
{
  if (filetypestr)
    {
      char *ftstr = filetypestr;

      // clang-format off
      if      (memcmp(filetypestr, "grb2",   4)  == 0) { ftstr += 4; DefaultFileType = CDI_FILETYPE_GRB2;}
      else if (memcmp(filetypestr, "grb1",   4)  == 0) { ftstr += 4; DefaultFileType = CDI_FILETYPE_GRB; }
      else if (memcmp(filetypestr, "grb",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_GRB; }
      else if (memcmp(filetypestr, "nc2",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_NC2; }
      else if (memcmp(filetypestr, "nc4c",   4)  == 0) { ftstr += 4; DefaultFileType = CDI_FILETYPE_NC4C;}
      else if (memcmp(filetypestr, "nc4",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_NC4; }
      else if (memcmp(filetypestr, "nc5",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_NC5; }
      else if (memcmp(filetypestr, "nc1",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_NC;  }
      else if (memcmp(filetypestr, "nczarr", 6)  == 0) { ftstr += 6; DefaultFileType = CDI_FILETYPE_NCZARR;  }
      else if (memcmp(filetypestr, "nc",     2)  == 0) { ftstr += 2; DefaultFileType = CDI_FILETYPE_NC2; }
      else if (memcmp(filetypestr, "srv",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_SRV; }
      else if (memcmp(filetypestr, "ext",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_EXT; }
      else if (memcmp(filetypestr, "ieg",    3)  == 0) { ftstr += 3; DefaultFileType = CDI_FILETYPE_IEG; }
      else
	{
	  fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
	  fprintf(stderr, "Available filetypes: grb, grb2, nc, nc2, nc4, nc4c, nc5, nczarr, srv, ext and ieg\n");
	  exit(EXIT_FAILURE);
	}
      // clang-format on

      if (DefaultFileType != CDI_UNDEFID && *ftstr != 0)
        {
          if (*ftstr == '_')
            {
              ftstr++;

              setDefaultDataType(ftstr);
            }
          else
            {
              fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
              fprintf(stderr, "Use format[_nbits] with:\n");
              fprintf(stderr, "    format = grb, grb2, nc, nc2, nc4, nc4c, nc5, nczarr, srv, ext or ieg\n");
              fprintf(stderr, "    nbits  = 32/64 for grb2/nc/nc2/nc4/nc4c/nc5/nczarr/srv/ext/ieg; 1 - 24 for grb/grb2\n");
              exit(EXIT_FAILURE);
            }
        }
    }
}

static int
handle_error(int cdiErrno, const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);

  printf("\n");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");

  va_end(args);

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  return cdiErrno;
}

static void
define_compress(const char *arg)
{
  size_t len = strlen(arg);

  if (strncmp(arg, "szip", len) == 0)
    {
      comptype = CDI_COMPRESS_SZIP;
    }
  else if (strncmp(arg, "jpeg", len) == 0)
    {
      comptype = CDI_COMPRESS_JPEG;
    }
  else if (strncmp(arg, "zip", 3) == 0)
    {
      comptype = CDI_COMPRESS_ZIP;
      complevel = (len == 5 && arg[3] == '_' && isdigit(arg[4])) ? atoi(&arg[4]) : 1;
    }
  else if (strncmp(arg, "zstd", 4) == 0)
    {
      if (filterId == 0)
        {
          filterId = 32015;
          numParams = 1;
          params[0] = (len >= 6 && len <= 7 && arg[4] == '_' && isdigit(arg[5])) ? atoi(&arg[5]) : 1;
        }
      else
        {
          fprintf(stderr, "Filter already defined!\n");
        }
    }
  else
    fprintf(stderr, "%s compression unsupported!\n", arg);
}

static void
define_filter(const char *arg)
{
  size_t len = strlen(arg);

  if (len > 0 && isdigit(arg[0]))
    {
      filterId = atoi(arg);
      const char *next;
      const char *parg = arg;
      while ((next = strchr(parg, ',')))
        {
          parg = ++next;
          len = strlen(parg);
          if (len > 0 && isdigit(parg[0]))
            {
              if (numParams >= maxParams)
                {
                  fprintf(stderr, "Too many filter parameter (max=%d)!\n", maxParams);
                  return;
                }
              params[numParams] = atoi(parg);
              numParams++;
            }
        }
    }
  else
    {
      fprintf(stderr, "Filter id missing!\n");
    }
}

int
main(int argc, char *argv[])
{
  char *rTable = NULL;
  char *wTable = NULL;
  bool preScan = false;
  int Move = 0;
  int Record = 0;
  int Variable = 0;
  int Debug = 0;
  int Vardis = 0;
  int Version = 0;
  int Longinfo = 0;
  int Shortinfo = 0;
  int varID;
  int itableID = CDI_UNDEFID, otableID = CDI_UNDEFID;
  int Info = 1;
  int NoInfo = 0;
  int numWorkerIn = 0;
  int numWorkerOut = 0;
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char *queryNames[256];
  int numQueryNames = 0;
  size_t queryCellidx[256];
  int numQueryCellidx = 0;
  int queryLevidx[256];
  int numQueryLevidx = 0;
  int queryStepidx[256];
  int numQueryStepidx = 0;

  Progname = strrchr(argv[0], '/');
  if (Progname == 0)
    Progname = argv[0];
  else
    Progname++;

  // clang-format off
  int c;
  while ((c = getopt(argc, argv, "b:C:F:f:i:L:o:N:S:t:w:z:cdEhlMmnRrsTvVxXZ")) != EOF)
    {
      switch (c)
	{
	case 'b': setDefaultDataType(optarg); break;
	case 'C': queryCellidx[numQueryCellidx++] = atol(optarg); break;
	case 'd': Debug = 1; break;
        case 'E': cdiDefGlobal("ECCODES_GRIB1", true); break;
	case 'f': setDefaultFileType(optarg); break;
	case 'h': usage(); exit (0);
	case 'i': numWorkerIn = atoi(optarg);  break;
	case 'o': numWorkerOut = atoi(optarg); break;
	case 'L': queryLevidx[numQueryLevidx++] = atoi(optarg); break;
	case 'l': Longinfo = 1; break;
	case 'M': cdiDefGlobal("HAVE_MISSVAL", 1); break;
	case 'm': Move = 1; break;
	case 'N': queryNames[numQueryNames++] = strdup(optarg); break;
	case 'n': Info = 0;  NoInfo = 1;  break;
	case 'R': cdiDefGlobal("REGULARGRID", 1); break;
	case 'r': Record = 1; break;
	case 'X': Variable = 1; break;
	case 'S': queryStepidx[numQueryStepidx++] = atoi(optarg); break;
	case 's': Shortinfo = 1; break;
	case 'T': preScan = true;  break;
	case 't': rTable = optarg;  break;
	case 'v': Vardis = 1; break;
	case 'V': Version = 1; break;
	case 'w': wTable = optarg; break;
	case 'x': datamode = SP_MODE; break;
	case 'z': define_compress(optarg); break;
	case 'F': define_filter(optarg); break;
	}
    }
  // clang-format on

  char *fname1 = (optind < argc) ? argv[optind++] : NULL;
  char *fname2 = (optind < argc) ? argv[optind++] : NULL;
  if (optind < argc) fprintf(stderr, "optind: %d argc: %d\n", optind, argc);

  if (Debug || Version) version();

  if (Debug) cdiDebug(Debug);

  if (rTable)
    {
      itableID = tableInq(-1, 0, rTable);
      if (itableID != CDI_UNDEFID) cdiDefTableID(itableID);
      otableID = itableID;
    }

  if (fname1 == NULL && !(Debug || Version))
    {
      usage();
      exit(0);
    }

  // for (int i = 0; i < numQueryNames; ++i) printf("queryName[%d] = %s\n", i+1, queryNames[i]);

  if (fname1)
    {
      SizeType nmiss;
      size_t datasize = 0;
      size_t maxlev = 0;
      int streamID2 = CDI_UNDEFID;
      int filetype;
      int nrecs;
      int levelID;
      int nts = 0;
      int recID;
      int taxisID2 = CDI_UNDEFID;
      int vlistID2 = CDI_UNDEFID;

      bool useQuery = (numQueryNames > 0) || (numQueryCellidx > 0) || (numQueryLevidx > 0) || (numQueryStepidx > 0);
      CdiQuery *query = NULL;
      if (useQuery)
        {
          query = cdiQueryCreate();
          if (numQueryNames) cdiQuerySetNames(query, numQueryNames, queryNames);
          if (numQueryCellidx) cdiQuerySetCellidx(query, numQueryCellidx, queryCellidx);
          if (numQueryLevidx) cdiQuerySetLevidx(query, numQueryLevidx, queryLevidx);
          if (numQueryStepidx) cdiQuerySetStepidx(query, numQueryStepidx, queryStepidx);
          if (Debug) cdiQueryPrint(query);
        }

      int streamID1 = useQuery ? streamOpenReadQuery(fname1, query) : streamOpenRead(fname1);
      if (streamID1 < 0) return handle_error(streamID1, "Open failed on %s", fname1);

      if (useQuery) cdiQueryPrintEntriesNotFound(query);
      if (query) cdiQueryDelete(query);

      if (numWorkerIn > 0) streamDefNumWorker(streamID1, numWorkerIn);

      if (preScan) fprintf(stdout, "numSteps=%d\n", streamInqNumSteps(streamID1));

      int vlistID1 = streamInqVlist(streamID1);

      if (Longinfo) vlistPrint(vlistID1);

      int nvars = vlistNvars(vlistID1);
      int taxisID1 = vlistInqTaxis(vlistID1);
      int ntsteps = vlistNtsteps(vlistID1);

      if (Debug) fprintf(stderr, "nvars   = %d\nntsteps = %d\n", nvars, ntsteps);

      if (fname2)
        {
          vlistID2 = vlistDuplicate(vlistID1);
          taxisID2 = taxisDuplicate(taxisID1);
          vlistDefTaxis(vlistID2, taxisID2);
        }

      for (varID = 0; varID < nvars; varID++)
        {
          int gridID = vlistInqVarGrid(vlistID1, varID);
          size_t gridsize = gridInqSize(gridID);
          if (gridsize > datasize) datasize = gridsize;
          int zaxisID = vlistInqVarZaxis(vlistID1, varID);
          size_t nlev = zaxisInqSize(zaxisID);
          if (nlev > maxlev) maxlev = nlev;
          if (fname2)
            {
              if (DefaultDataType != CDI_UNDEFID) vlistDefVarDatatype(vlistID2, varID, DefaultDataType);
            }
        }

      if (fname2)
        {
          Info = 0;
          filetype = streamInqFiletype(streamID1);

          if (DefaultFileType != CDI_UNDEFID) filetype = DefaultFileType;

          streamID2 = streamOpenWrite(fname2, filetype);
          if (streamID2 < 0) return handle_error(streamID2, "Open failed on %s", fname2);

          int maxSteps = vlistNtsteps(vlistID1);
          if (filetype == CDI_FILETYPE_NCZARR && maxSteps >= 0) streamDefMaxSteps(streamID2, maxSteps);
          if (numWorkerOut > 0) streamDefNumWorker(streamID2, numWorkerOut);

          if (DefaultByteorder != CDI_UNDEFID) streamDefByteorder(streamID2, DefaultByteorder);

          if (filterId != 0)
            {
              streamDefFilter(streamID2, filterId, numParams, params);
            }

          if (comptype != CDI_COMPRESS_NONE)
            {
              streamDefCompType(streamID2, comptype);
              streamDefCompLevel(streamID2, complevel);
            }

          streamDefVlist(streamID2, vlistID2);

          if (otableID == CDI_UNDEFID) otableID = itableID;
        }

      if (vlistNumber(vlistID1) != CDI_REAL) datasize *= 2;
      if (Variable) datasize *= maxlev;
      double *data = (double *) malloc(datasize * sizeof(double));

      // nts = cdiInqTimeSize(streamID1);
      if (Debug) printf("nts = %d streamID1 = %d, streamID2 = %d\n", nts, streamID1, streamID2);

      if (Shortinfo)
        {
          Info = 0;
          print_short_info(streamID1, vlistID1, Vardis);
        }

      int tsID = 0;
      if (Info || fname2 || NoInfo)
        while ((nrecs = streamInqTimestep(streamID1, tsID)) > 0)
          {
            if (fname2 /* && ntsteps != 0*/)
              {
                taxisCopyTimestep(taxisID2, taxisID1);
                streamDefTimestep(streamID2, tsID);
              }

            CdiDateTime vdatetime = taxisInqVdatetime(taxisID1);

            if (Record)
              {
                for (recID = 0; recID < nrecs; recID++)
                  {
                    streamInqRecord(streamID1, &varID, &levelID);
                    streamReadRecord(streamID1, data, &nmiss);

                    int number = vlistInqVarNumber(vlistID1, varID);
                    int gridID = vlistInqVarGrid(vlistID1, varID);
                    int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                    int param = vlistInqVarParam(vlistID1, varID);

                    cdiParamToString(param, paramstr, sizeof(paramstr));

                    if (Vardis)
                      vlistInqVarName(vlistID1, varID, varname);
                    else
                      strcpy(varname, paramstr);

                    // printf("varID=%d, param=%d, gridID=%d, zaxisID=%d levelID=%d\n", varID, param, gridID, zaxisID, levelID);

                    SizeType gridsize = gridInqSize(gridID);
                    double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) : levelID + 1;
                    double missval = vlistInqVarMissval(vlistID1, varID);

                    if (Info)
                      printInfo(vdatetime, varname, level, (size_t) gridsize, number, (size_t) nmiss, missval, data, Vardis);

                    if (fname2)
                      {
                        streamDefRecord(streamID2, varID, levelID);
                        if (Move)
                          streamCopyRecord(streamID2, streamID1);
                        else
                          streamWriteRecord(streamID2, data, nmiss);
                      }
                  }
              }
            else if (Variable)
              {
                for (varID = 0; varID < nvars; varID++)
                  {
                    if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT && tsID > 0) continue;

                    int number = vlistInqVarNumber(vlistID1, varID);
                    int gridID = vlistInqVarGrid(vlistID1, varID);
                    int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                    int param = vlistInqVarParam(vlistID1, varID);

                    cdiParamToString(param, paramstr, sizeof(paramstr));

                    if (Vardis)
                      vlistInqVarName(vlistID1, varID, varname);
                    else
                      strcpy(varname, paramstr);

                    if (Debug) fprintf(stdout, "varID = %d param = %d gridID = %d zaxisID = %d\n", varID, param, gridID, zaxisID);

                    size_t gridsize = gridInqSize(gridID);
                    double missval = vlistInqVarMissval(vlistID1, varID);

                    streamReadVar(streamID1, varID, data, &nmiss);

                    int nlevs = zaxisInqSize(zaxisID);
                    for (levelID = 0; levelID < nlevs; levelID++)
                      {
                        double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) : levelID + 1;
                        size_t offset = levelID * gridsize;
                        if (Info) printInfo(vdatetime, varname, level, gridsize, number, nmiss, missval, data + offset, Vardis);
                      }

                    if (fname2) streamWriteVar(streamID2, varID, data, nmiss);
                  }
              }
            else
              {
                for (varID = 0; varID < nvars; varID++)
                  {
                    if (vlistInqVarTimetype(vlistID1, varID) == TIME_CONSTANT && tsID > 0) continue;

                    int number = vlistInqVarNumber(vlistID1, varID);
                    int gridID = vlistInqVarGrid(vlistID1, varID);
                    int zaxisID = vlistInqVarZaxis(vlistID1, varID);
                    int param = vlistInqVarParam(vlistID1, varID);

                    cdiParamToString(param, paramstr, sizeof(paramstr));

                    if (Vardis)
                      vlistInqVarName(vlistID1, varID, varname);
                    else
                      strcpy(varname, paramstr);

                    if (Debug) fprintf(stdout, "varID = %d param = %d gridID = %d zaxisID = %d\n", varID, param, gridID, zaxisID);

                    size_t gridsize = gridInqSize(gridID);
                    double missval = vlistInqVarMissval(vlistID1, varID);

                    int nlevs = zaxisInqSize(zaxisID);
                    for (levelID = 0; levelID < nlevs; levelID++)
                      {
                        double level = zaxisInqLevels(zaxisID, NULL) ? zaxisInqLevel(zaxisID, levelID) : levelID + 1;
                        streamReadVarSlice(streamID1, varID, levelID, data, &nmiss);

                        if (Info) printInfo(vdatetime, varname, level, gridsize, number, nmiss, missval, data, Vardis);

                        if (fname2) streamWriteVarSlice(streamID2, varID, levelID, data, nmiss);
                      }
                  }
              }

            tsID++;
          }

      free(data);

      // fprintf(stderr, "%ld\n", (long) streamNvals(streamID1));

      if (fname2)
        {
          streamClose(streamID2);
          vlistDestroy(vlistID2);
          taxisDestroy(taxisID2);
        }
      streamClose(streamID1);
    }

  for (int i = 0; i < numQueryNames; ++i) free(queryNames[i]);

  if (wTable) tableWrite(wTable, itableID);

  return 0;
}
