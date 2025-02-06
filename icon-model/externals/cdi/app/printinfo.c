#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include <cdi.h>
#include "cdi_uuid.h"
#include "cdi_int.h"

#include "printinfo.h"

#define DATE_FORMAT "%5.4d-%2.2d-%2.2d"
#define TIME_FORMAT "%2.2d:%2.2d:%2.2d"

void
datetime2str(CdiDateTime dt, char *datetimestr, int maxlen)
{
  snprintf(datetimestr, maxlen, DATE_FORMAT "T" TIME_FORMAT, dt.date.year, dt.date.month, dt.date.day, dt.time.hour, dt.time.minute,
           dt.time.second);
}

void
date2str(CdiDate date, char *datestr, int maxlen)
{
  snprintf(datestr, maxlen, DATE_FORMAT, date.year, date.month, date.day);
}

void
time2str(CdiTime time, char *timestr, int maxlen)
{
  static bool readEnv = true;
  static int msDigitsNum = 0;
  if (readEnv)
    {
      readEnv = false;
      char *envString = getenv("CDI_MS_DIGITS");
      if (envString)
        {
          int ival = atol(envString);
          if (ival > 0) msDigitsNum = ival;
          if (ival > 3) msDigitsNum = 3;
        }
    }

  if (msDigitsNum)
    snprintf(timestr, maxlen, "%2.2d:%2.2d:%0*.*f", time.hour, time.minute, msDigitsNum + 3, msDigitsNum,
             time.second + time.ms / 1000.0);
  else
    snprintf(timestr, maxlen, TIME_FORMAT, time.hour, time.minute, time.second);
}

const char *
comptype_to_name(int compType)
{
  switch (compType)
    {
    case CDI_COMPRESS_SZIP: return "szip";
    case CDI_COMPRESS_AEC: return "aec";
    case CDI_COMPRESS_ZIP: return "zip";
    case CDI_COMPRESS_JPEG: return "jpeg";
    case CDI_COMPRESS_FILTER: return "filter";
    }
  return " ";
}

void
printFiletype(int streamID, int vlistID)
{
  int filetype = streamInqFiletype(streamID);

  // clang-format off
  switch (filetype)
    {
    case CDI_FILETYPE_GRB:    printf("GRIB");  break;
    case CDI_FILETYPE_GRB2:   printf("GRIB2");  break;
    case CDI_FILETYPE_NC:     printf("NetCDF");  break;
    case CDI_FILETYPE_NC2:    printf("NetCDF2");  break;
    case CDI_FILETYPE_NC4:    printf("NetCDF4");  break;
    case CDI_FILETYPE_NC4C:   printf("NetCDF4 classic");  break;
    case CDI_FILETYPE_NC5:    printf("NetCDF5");  break;
    case CDI_FILETYPE_NCZARR: printf("NCZarr");  break;
    case CDI_FILETYPE_SRV:    printf("SERVICE");  break;
    case CDI_FILETYPE_EXT:    printf("EXTRA");  break;
    case CDI_FILETYPE_IEG:    printf("IEG");  break;
    default: printf("  File format: unsupported filetype %d" , filetype);  break;
    }

  if (filetype == CDI_FILETYPE_SRV || filetype == CDI_FILETYPE_EXT || filetype == CDI_FILETYPE_IEG)
    {
      switch (streamInqByteorder(streamID))
	{
	case CDI_BIGENDIAN:    printf("  BIGENDIAN");  break;
	case CDI_LITTLEENDIAN: printf("  LITTLEENDIAN");  break;
	default: printf("  byteorder: %d undefined", streamInqByteorder(streamID));  break;
	}
    }
  // clang-format on

  int nvars = vlistNvars(vlistID);
  const int comps[] = { CDI_COMPRESS_SZIP, CDI_COMPRESS_AEC, CDI_COMPRESS_ZIP, CDI_COMPRESS_JPEG, CDI_COMPRESS_FILTER };
  unsigned kk = 0;
  for (unsigned k = 0; k < sizeof(comps) / sizeof(int); ++k)
    for (int varID = 0; varID < nvars; varID++)
      {
        int comptype = vlistInqVarCompType(vlistID, varID);
        if (comptype == comps[k])
          {
            printf("%c%s", (kk++ == 0) ? ' ' : '/', comptype_to_name(comptype));
            break;
          }
      }

  printf("\n");
}

static void
print_xvals(int gridID, int dig)
{
  size_t xsize = gridInqXsize(gridID);
  if (xsize > 0 && gridInqXvals(gridID, NULL))
    {
      char xname[CDI_MAX_NAME], xunits[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, xname, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, xunits, &length);

      const double xfirst = gridInqXval(gridID, 0);
      const double xlast = gridInqXval(gridID, xsize - 1);
      const double xinc = gridInqXinc(gridID);
      fprintf(stdout, "%33s : %.*g", xname, dig, xfirst);
      if (xsize > 1)
        {
          fprintf(stdout, " to %.*g", dig, xlast);
          if (IS_NOT_EQUAL(xinc, 0)) fprintf(stdout, " by %.*g", dig, xinc);
        }
      fprintf(stdout, " %s%s\n", xunits, gridIsCircular(gridID) ? "  circular" : "");
    }
}

static void
print_yvals(int gridID, int dig)
{
  size_t ysize = gridInqYsize(gridID);
  if (ysize > 0 && gridInqYvals(gridID, NULL))
    {
      char yname[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, yname, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, yunits, &length);

      const double yfirst = gridInqYval(gridID, 0);
      const double ylast = gridInqYval(gridID, ysize - 1);
      const double yinc = gridInqYinc(gridID);
      fprintf(stdout, "%33s : %.*g", yname, dig, yfirst);
      if (ysize > 1)
        {
          int gridtype = gridInqType(gridID);
          fprintf(stdout, " to %.*g", dig, ylast);
          if (IS_NOT_EQUAL(yinc, 0) && gridtype != GRID_GAUSSIAN && gridtype != GRID_GAUSSIAN_REDUCED)
            fprintf(stdout, " by %.*g", dig, yinc);
        }
      fprintf(stdout, " %s\n", yunits);
    }
}

static void
print_xyvals2D(int gridID, int dig)
{
  if (gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL))
    {
      char xname[CDI_MAX_NAME], yname[CDI_MAX_NAME], xunits[CDI_MAX_NAME], yunits[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, xname, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, yname, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, xunits, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, yunits, &length);

      size_t gridsize = gridInqSize(gridID);
      double *xvals2D = (double *) malloc(gridsize * sizeof(double));
      double *yvals2D = (double *) malloc(gridsize * sizeof(double));

      gridInqXvals(gridID, xvals2D);
      gridInqYvals(gridID, yvals2D);

      double xfirst = xvals2D[0];
      double xlast = xvals2D[0];
      double yfirst = yvals2D[0];
      double ylast = yvals2D[0];
      for (size_t i = 1; i < gridsize; i++)
        {
          if (xvals2D[i] < xfirst) xfirst = xvals2D[i];
          if (xvals2D[i] > xlast) xlast = xvals2D[i];
          if (yvals2D[i] < yfirst) yfirst = yvals2D[i];
          if (yvals2D[i] > ylast) ylast = yvals2D[i];
        }

      double xinc = 0;
      double yinc = 0;
      int gridtype = gridInqType(gridID);
      if (gridtype == GRID_CURVILINEAR)
        {
          size_t xsize = gridInqXsize(gridID);
          size_t ysize = gridInqYsize(gridID);
          if (xsize > 1)
            {
              double *xvals = (double *) malloc((size_t) xsize * sizeof(double));
              for (size_t i = 0; i < xsize; ++i) xvals[i] = xvals2D[i];
              xinc = fabs(xvals[xsize - 1] - xvals[0]) / (xsize - 1);
              for (size_t i = 1; i < xsize; i++)
                if (fabs(fabs(xvals[i - 1] - xvals[i]) - xinc) > 0.01 * xinc)
                  {
                    xinc = 0;
                    break;
                  }
              free(xvals);
              if (IS_NOT_EQUAL(xinc, 0))
                {
                  for (size_t i = 1; i < ysize; i++)
                    if (IS_NOT_EQUAL(xvals2D[i * xsize], xvals2D[0])
                        || IS_NOT_EQUAL(xvals2D[(i + 1) * xsize - 1], xvals2D[xsize - 1]))
                      {
                        xinc = 0;
                        break;
                      }
                }
            }
          if (ysize > 1)
            {
              double *yvals = (double *) malloc((size_t) ysize * sizeof(double));
              for (size_t i = 0; i < ysize; ++i) yvals[i] = yvals2D[i * xsize];
              yinc = fabs(yvals[ysize - 1] - yvals[0]) / (ysize - 1);
              for (size_t i = 1; i < ysize; i++)
                if (fabs(fabs(yvals[i - 1] - yvals[i]) - yinc) > 0.01 * yinc)
                  {
                    yinc = 0;
                    break;
                  }
              free(yvals);
              if (IS_NOT_EQUAL(yinc, 0))
                {
                  for (size_t i = 1; i < xsize; i++)
                    if (IS_NOT_EQUAL(yvals2D[i], yvals2D[0])
                        || IS_NOT_EQUAL(yvals2D[(ysize - 1) * xsize + i], yvals2D[(ysize - 1) * xsize]))
                      {
                        yinc = 0;
                        break;
                      }
                }
            }
        }

      fprintf(stdout, "%33s : %.*g", xname, dig, xfirst);
      if (gridsize > 1) fprintf(stdout, " to %.*g", dig, xlast);
      if (IS_NOT_EQUAL(xinc, 0)) fprintf(stdout, " by %.*g", dig, xinc);
      fprintf(stdout, " %s%s\n", xunits, gridIsCircular(gridID) ? "  circular" : "");
      fprintf(stdout, "%33s : %.*g", yname, dig, yfirst);
      if (gridsize > 1) fprintf(stdout, " to %.*g", dig, ylast);
      if (IS_NOT_EQUAL(yinc, 0)) fprintf(stdout, " by %.*g", dig, yinc);
      fprintf(stdout, " %s\n", yunits);

      free(xvals2D);
      free(yvals2D);
    }
}

static void
printGridNp(int gridtype, int gridID, size_t gridsize, size_t xsize, size_t ysize)
{
  fprintf(stdout, "points=%zu", gridsize);
  if (gridtype == GRID_GAUSSIAN_REDUCED)
    fprintf(stdout, "  nlat=%zu", ysize);
  else if (xsize && ysize)
    fprintf(stdout, " (%zux%zu)", xsize, ysize);

  int numLPE = gridInqNP(gridID);
  if (numLPE > 0)
    {
      if (gridtype == GRID_GAUSSIAN) fprintf(stdout, "  F%d", numLPE);
      if (gridtype == GRID_GAUSSIAN_REDUCED) fprintf(stdout, "  N%d", numLPE);
    }

  fprintf(stdout, "\n");
}

static void
printGridInfoKernel(int gridID, int index, bool lproj)
{
  int gridtype = gridInqType(gridID);

  if (lproj && gridtype != GRID_PROJECTION)
    fprintf(stderr, "Internal problem (%s): sub grid not equal GRID_PROJECTION!\n", __func__);

  int trunc = gridInqTrunc(gridID);
  size_t gridsize = gridInqSize(gridID);
  size_t xsize = gridInqXsize(gridID);
  size_t ysize = gridInqYsize(gridID);

  // int prec     = gridInqDatatype(gridID);
  // int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;
  int dig = 7;

  if (!lproj)
    {
      fprintf(stdout, "  %4d : ", index + 1);
      fprintf(stdout, "%-24s", gridNamePtr(gridtype));
      fprintf(stdout, " : ");
    }

  if (gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION || gridtype == GRID_GENERIC || gridtype == GRID_CHARXY
      || gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      if (!lproj)
        {
          printGridNp(gridtype, gridID, gridsize, xsize, ysize);
        }

      char name[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, name, &length);
      if (gridtype == GRID_PROJECTION || name[0])
        {
          if (name[0] == 0) strcpy(name, "undefined");
          fprintf(stdout, "         %24s", "mapping");
          fprintf(stdout, " : ");
          fprintf(stdout, "%s\n", name);
        }

      print_xvals(gridID, dig);
      print_yvals(gridID, dig);

      if (gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL))
        {
          fprintf(stdout, "%33s :%s%s%s\n", "available",
                  (gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL) ? " cellbounds" : ""),
                  gridHasArea(gridID) ? " area" : "", gridInqMask(gridID, NULL) ? " mask" : "");
        }
    }
  else if (gridtype == GRID_SPECTRAL)
    {
      fprintf(stdout, "points=%zu  nsp=%zu  T%d%s\n", gridsize, gridsize / 2, trunc,
              gridInqComplexPacking(gridID) ? "  complexPacking" : "");
    }
  else if (gridtype == GRID_FOURIER)
    {
      fprintf(stdout, "points=%zu  nfc=%zu  T%d\n", gridsize, gridsize / 2, trunc);
    }
  else if (gridtype == GRID_GME)
    {
      int nd, ni, ni2, ni3;
      gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
      fprintf(stdout, "points=%zu  nd=%d  ni=%d\n", gridsize, nd, ni);
    }
  else if (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
    {
      if (gridtype == GRID_CURVILINEAR)
        fprintf(stdout, "points=%zu (%zux%zu)", gridsize, xsize, ysize);
      else
        fprintf(stdout, "points=%zu", gridsize);

      if (gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) > 0) fprintf(stdout, "  nvertex=%d", gridInqNvertex(gridID));

      fprintf(stdout, "\n");

      if (gridtype == GRID_UNSTRUCTURED)
        {
          int number = gridInqNumber(gridID);
          int position = gridInqPosition(gridID);
          if (number > 0) fprintf(stdout, "%33s : number=%d  position=%d\n", "grid", number, position);

          if (gridInqReference(gridID, NULL))
            {
              char reference_link[8192];
              gridInqReference(gridID, reference_link);
              fprintf(stdout, "%33s : %s\n", "uri", reference_link);
            }
        }

      print_xyvals2D(gridID, dig);
    }
  else /* if ( gridtype == GRID_GENERIC ) */
    {
      if (ysize == 0)
        fprintf(stdout, "points=%zu\n", gridsize);
      else
        fprintf(stdout, "points=%zu (%zux%zu)\n", gridsize, xsize, ysize);
    }

  if (gridtype == GRID_CURVILINEAR || gridtype == GRID_UNSTRUCTURED)
    {
      if (gridHasArea(gridID) || gridInqXbounds(gridID, NULL) || gridInqYbounds(gridID, NULL))
        {
          fprintf(stdout, "%33s :%s%s%s\n", "available",
                  (gridInqXbounds(gridID, NULL) && gridInqYbounds(gridID, NULL)) ? " cellbounds" : "",
                  gridHasArea(gridID) ? " area" : "", gridInqMask(gridID, NULL) ? " mask" : "");
        }
    }

  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  gridInqUUID(gridID, uuidOfHGrid);
  if (!cdiUUIDIsNull(uuidOfHGrid))
    {
      char uuidOfHGridStr[37];
      cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
      if (uuidOfHGridStr[0] != 0 && strlen(uuidOfHGridStr) == 36)
        {
          fprintf(stdout, "%33s : %s\n", "uuid", uuidOfHGridStr);
        }
    }
}

void
printGridInfo(int vlistID)
{
  int ngrids = vlistNgrids(vlistID);
  for (int index = 0; index < ngrids; index++)
    {
      int gridID = vlistGrid(vlistID, index);
      printGridInfoKernel(gridID, index, false);
      int projID = gridInqProj(gridID);
      if (projID != CDI_UNDEFID) printGridInfoKernel(projID, index, true);
    }
}

static void
printZaxisBoundsInfo(int zaxisID, int dig, int levelsize, const double zinc, const char *zunits)
{
  double level1 = zaxisInqLbound(zaxisID, 0);
  double level2 = zaxisInqUbound(zaxisID, 0);
  fprintf(stdout, "%33s : %.*g-%.*g", "bounds", dig, level1, dig, level2);
  if (levelsize > 1)
    {
      level1 = zaxisInqLbound(zaxisID, levelsize - 1);
      level2 = zaxisInqUbound(zaxisID, levelsize - 1);
      fprintf(stdout, " to %.*g-%.*g", dig, level1, dig, level2);
      if (IS_NOT_EQUAL(zinc, 0)) fprintf(stdout, " by %.*g", dig, zinc);
    }
  fprintf(stdout, " %s\n", zunits);
}

static void
printZaxisLevelInfo(int levelsize, int zaxisID, int zaxistype, double zinc, int dig, const char *zname, const char *zunits)
{
  double *levels = (double *) malloc((size_t) levelsize * sizeof(double));
  zaxisInqLevels(zaxisID, levels);

  if (!(zaxistype == ZAXIS_SURFACE && levelsize == 1 && fabs(levels[0]) <= 0))
    {
      const double zfirst = levels[0];
      const double zlast = levels[levelsize - 1];
      if (levelsize > 2)
        {
          zinc = (levels[levelsize - 1] - levels[0]) / (levelsize - 1);
          for (int levelID = 2; levelID < levelsize; ++levelID)
            if (fabs(fabs(levels[levelID] - levels[levelID - 1]) - zinc) > 0.001 * zinc)
              {
                zinc = 0;
                break;
              }
        }

      fprintf(stdout, "%33s : %.*g", zname, dig, zfirst);
      if (levelsize > 1)
        {
          fprintf(stdout, " to %.*g", dig, zlast);
          if (IS_NOT_EQUAL(zinc, 0)) fprintf(stdout, " by %.*g", dig, zinc);
        }
      fprintf(stdout, " %s\n", zunits);
    }

  free(levels);
}

static void
printZaxisHybridInfo(int zaxisID)
{
  char psname[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_PSNAME, psname, &length);
  int vctsize = zaxisInqVctSize(zaxisID);
  if (vctsize || psname[0])
    {
      fprintf(stdout, "%33s :%s%s%s\n", "available", vctsize ? " vct" : "", psname[0] ? "  ps: " : "", psname);
    }
}

static void
printZaxisGenericInfo(int ltype, int zaxistype, const char *zaxisname)
{
  if (zaxistype == ZAXIS_GENERIC && ltype != 0)
    fprintf(stdout, "%-12s (ltype=%3d)", zaxisname, ltype);
  else
    fprintf(stdout, "%-24s", zaxisname);
}

static void
printZaxisReferenceInfo(int zaxisID)
{
  int referenceNumber = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_NUMBEROFVGRIDUSED, &referenceNumber);
  if (referenceNumber > 0)
    {
      fprintf(stdout, "%33s : number=%d\n", "zaxis", referenceNumber);
    }

  unsigned char uuidOfVGrid[CDI_UUID_SIZE];
  int length = CDI_UUID_SIZE;
  cdiInqKeyBytes(zaxisID, CDI_GLOBAL, CDI_KEY_UUID, uuidOfVGrid, &length);
  if (!cdiUUIDIsNull(uuidOfVGrid))
    {
      char uuidOfVGridStr[37];
      cdiUUID2Str(uuidOfVGrid, uuidOfVGridStr);
      if (uuidOfVGridStr[0] != 0 && strlen(uuidOfVGridStr) == 36)
        {
          fprintf(stdout, "%33s : %s\n", "uuid", uuidOfVGridStr);
        }
    }
}

void
printZaxisInfo(int vlistID)
{
  char zaxisname[CDI_MAX_NAME], zname[CDI_MAX_NAME], zunits[CDI_MAX_NAME];

  int nzaxis = vlistNzaxis(vlistID);
  for (int index = 0; index < nzaxis; index++)
    {
      int dig = 7;
      //______________________________________________________--
      double zinc = 0;
      int zaxisID = vlistZaxis(vlistID, index);
      int zaxistype = zaxisInqType(zaxisID);
      int ltype = 0;
      cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);
      int levelsize = zaxisInqSize(zaxisID);
      // int prec      = zaxisInqDatatype(zaxisID);
      // int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;

      zaxisName(zaxistype, zaxisname);
      int length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, zname, &length);
      length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, zunits, &length);
      zunits[12] = 0;

      fprintf(stdout, "  %4d : ", vlistZaxisIndex(vlistID, zaxisID) + 1);

      printZaxisGenericInfo(ltype, zaxistype, zaxisname);

      fprintf(stdout, " :");

      const bool zscalar = (levelsize == 1) ? zaxisInqScalar(zaxisID) : false;
      fprintf(stdout, " levels=%d%s\n", levelsize, zscalar ? "  scalar" : "");

      if (zaxisInqLevels(zaxisID, NULL))
        {
          printZaxisLevelInfo(levelsize, zaxisID, zaxistype, zinc, dig, zname, zunits);
        }

      if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
        {
          printZaxisBoundsInfo(zaxisID, dig, levelsize, zinc, zunits);
        }

      if (zaxistype == ZAXIS_HYBRID) printZaxisHybridInfo(zaxisID);

      if (zaxistype == ZAXIS_REFERENCE) printZaxisReferenceInfo(zaxisID);
    }
}

void
printSubtypeInfo(int vlistID)
{
  int nsubtypes = vlistNsubtypes(vlistID);
  for (int index = 0; index < nsubtypes; index++)
    {
      int subtypeID = vlistSubtype(vlistID, index);
      int subtypesize = subtypeInqSize(subtypeID);
      // subtypePrint(subtypeID);
      fprintf(stdout, "  %4d : %-24s : ntiles=%d\n", vlistSubtypeIndex(vlistID, subtypeID) + 1, "tiles", subtypesize);
    }
}

static int
printDateTime(int ntimeout, CdiDateTime vdatetime)
{
  if (ntimeout == 4)
    {
      ntimeout = 0;
      fprintf(stdout, "\n");
    }

  char vdatestr[32], vtimestr[32];
  date2str(vdatetime.date, vdatestr, sizeof(vdatestr));
  time2str(vdatetime.time, vtimestr, sizeof(vtimestr));

  fprintf(stdout, " %s %s", vdatestr, vtimestr);

  return ++ntimeout;
}

#define NUM_TIMESTEP 60
#define MAX_DOTS 80

static int
printDot(int ndotout, int *nfact, int *ncout)
{
  // printf("ncout %d %d %d\n",*ncout, (*ncout)%(*nfact), *nfact);
  if ((*ncout) % (*nfact) == 0)
    {
      if (ndotout == MAX_DOTS)
        {
          *ncout = 0;
          ndotout = 0;
          fprintf(stdout, "\n   ");
          (*nfact) *= 10;
        }

      fprintf(stdout, ".");
      fflush(stdout);
      ndotout++;
    }

  (*ncout)++;

  return ndotout;
}

void
printTimesteps(int streamID, int taxisID, int verbose)
{
  struct DateTimeEntry
  {
    CdiDateTime dt;
    struct DateTimeEntry *next;
  };
  struct DateTimeEntry dateTimeList[NUM_TIMESTEP];
  struct DateTimeEntry *dateTimeNext = dateTimeList;

  for (int i = 0; i < NUM_TIMESTEP - 1; ++i) dateTimeList[i].next = &dateTimeList[i + 1];
  dateTimeList[NUM_TIMESTEP - 1].next = &dateTimeList[0];

  int ntimeout = 0;
  int ndotout = 0;
  int nvdatetime = 0;
  int ncout = 0;
  int nfact = 1;
  int tsID = 0;

  while (true)
    {
      int nrecs = streamInqTimestep(streamID, tsID);
      if (nrecs == 0) break;

      const CdiDateTime vdatetime = taxisInqVdatetime(taxisID);

      if (verbose || tsID < NUM_TIMESTEP)
        {
          ntimeout = printDateTime(ntimeout, vdatetime);
        }
      else
        {
          if (tsID == 2 * NUM_TIMESTEP) fprintf(stdout, "\n   ");
          if (tsID >= 2 * NUM_TIMESTEP) ndotout = printDot(ndotout, &nfact, &ncout);

          if (nvdatetime < NUM_TIMESTEP)
            {
              dateTimeList[nvdatetime].dt = vdatetime;
              nvdatetime++;
            }
          else
            {
              dateTimeNext->dt = vdatetime;
              dateTimeNext = dateTimeNext->next;
            }
        }

      tsID++;
    }

  if (nvdatetime)
    {
      fprintf(stdout, "\n");

      ntimeout = 0;
      int toff = 0;
      if (tsID > 2 * NUM_TIMESTEP)
        {
          toff = tsID % 4;
          if (toff > 0) toff = 4 - toff;
          for (int i = 0; i < toff; ++i) dateTimeNext = dateTimeNext->next;
        }
      for (int i = toff; i < nvdatetime; ++i)
        {
          ntimeout = printDateTime(ntimeout, dateTimeNext->dt);
          dateTimeNext = dateTimeNext->next;
        }
    }
}
