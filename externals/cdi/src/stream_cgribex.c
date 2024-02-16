#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdio.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "file.h"
#include "varscan.h"
#include "stream_scan.h"
#include "stream_grb.h"
#include "stream_cgribex.h"

#ifdef HAVE_LIBCGRIBEX

#include "cgribex.h"

typedef struct
{
  int sec0[2];
  int sec1[1024];
  size_t sec2len;
  int *sec2;
  int sec3[2];
  int sec4[512];
  double fsec2[512];
  double fsec3[2];
} cgribexrec_t;

typedef struct
{
  int param;
  int level1;
  int level2;
  int ltype;
  int tsteptype;
} compvar_t;

typedef struct
{
  void *gribbuffer;
  size_t gribbuffersize;
  unsigned char *pds;
  unsigned char *gds;
  unsigned char *bms;
  unsigned char *bds;
} cgribex_handle;

static void
fill_intarr(int *iarr, int val, int n)
{
  for (int i = 0; i < n; ++i) iarr[i] = val;
}

static void
cgribexInit(cgribexrec_t *cgribexp)
{
  cgribexp->sec2len = 4096;
  cgribexp->sec2 = (int *) Malloc(cgribexp->sec2len * sizeof(int));
}

void *
cgribexNew(void)
{
  cgribexrec_t *cgribexp = (cgribexrec_t *) Malloc(sizeof(cgribexrec_t));
  cgribexInit(cgribexp);
  return (void *) cgribexp;
}

void
cgribexDelete(void *cgribex)
{
  cgribexrec_t *cgribexp = (cgribexrec_t *) cgribex;
  if (cgribexp)
    {
      if (cgribexp->sec2) Free(cgribexp->sec2);
      Free(cgribexp);
    }
}

#ifdef __cplusplus
extern "C"
{
#endif
  int grib1Sections(unsigned char *gribbuffer, long gribbufsize, unsigned char **pdsp, unsigned char **gdsp, unsigned char **bmsp,
                    unsigned char **bdsp, long *gribrecsize);
#ifdef __cplusplus
}
#endif

static size_t
cgribexSection2Length(void *gribbuffer, size_t gribbuffersize)
{
  long sec2len = 0;

  if (gribbuffersize && gribbuffer)
    {
      unsigned char *pds = NULL, *gds = NULL, *bms = NULL, *bds = NULL;
      long gribrecsize;
      int status = grib1Sections((unsigned char *) gribbuffer, (long) gribbuffersize, &pds, &gds, &bms, &bds, &gribrecsize);
      if (status >= 0 && gds != NULL) sec2len = (unsigned) ((gds[0] << 16) + (gds[1] << 8) + (gds[3]));
    }

  return sec2len;
}

void *
cgribex_handle_new_from_meassage(void *gribbuffer, size_t gribbuffersize)
{
  cgribex_handle *gh = (cgribex_handle *) Malloc(sizeof(cgribex_handle));
  gh->gribbuffer = NULL;
  gh->gribbuffersize = 0;
  gh->pds = NULL;

  if (gribbuffersize && gribbuffer)
    {
      unsigned char *pds = NULL, *gds = NULL, *bms = NULL, *bds = NULL;
      long gribrecsize;
      int status = grib1Sections((unsigned char *) gribbuffer, (long) gribbuffersize, &pds, &gds, &bms, &bds, &gribrecsize);
      if (status >= 0)
        {
          gh->gribbuffer = gribbuffer;
          gh->gribbuffersize = gribbuffersize;
          gh->pds = pds;
          gh->gds = gds;
          gh->bms = bms;
          gh->bds = bds;
        }
    }

  return (void *) gh;
}

void
cgribex_handle_delete(void *gh)
{
  if (gh) Free(gh);
}

static int
cgribexGetGridType(int *isec2)
{
  int gridtype = GRID_GENERIC;

  // clang-format off
  switch (ISEC2_GridType)
    {
    case  GRIB1_GTYPE_LATLON:     { gridtype = GRID_LONLAT;     break; }
    case  GRIB1_GTYPE_LATLON_ROT: { gridtype = GRID_PROJECTION; break; }
    case  GRIB1_GTYPE_LCC:        { gridtype = CDI_PROJ_LCC;    break; }
    case  GRIB1_GTYPE_GAUSSIAN:   { gridtype = ISEC2_Reduced ? GRID_GAUSSIAN_REDUCED : GRID_GAUSSIAN; break; }
    case  GRIB1_GTYPE_SPECTRAL:   { gridtype = GRID_SPECTRAL;   break; }
    case  GRIB1_GTYPE_GME:        { gridtype = GRID_GME;        break; }
    }
  // clang-format on

  return gridtype;
}

static bool
cgribexGetIsRotated(int *isec2)
{
  return (ISEC2_GridType == GRIB1_GTYPE_LATLON_ROT);
}

static bool
cgribexGetZaxisHasBounds(int grb_ltype)
{
  // clang-format off
  switch (grb_ltype)
    {
    case GRIB1_LTYPE_SIGMA_LAYER:
    case GRIB1_LTYPE_HYBRID_LAYER:
    case GRIB1_LTYPE_LANDDEPTH_LAYER: return true;
    }
  // clang-format on

  return false;
}

static int
cgribexGetTimeUnit(int *isec1)
{
  int timeunit = TUNIT_HOUR;
  static bool lprint = true;

  // clang-format off
  switch ( ISEC1_TimeUnit )
    {
    case ISEC1_TABLE4_MINUTE:    timeunit = TUNIT_MINUTE;    break;
    case ISEC1_TABLE4_QUARTER:   timeunit = TUNIT_QUARTER;   break;
    case ISEC1_TABLE4_30MINUTES: timeunit = TUNIT_30MINUTES; break;
    case ISEC1_TABLE4_HOUR:      timeunit = TUNIT_HOUR;      break;
    case ISEC1_TABLE4_3HOURS:    timeunit = TUNIT_3HOURS;    break;
    case ISEC1_TABLE4_6HOURS:    timeunit = TUNIT_6HOURS;    break;
    case ISEC1_TABLE4_12HOURS:   timeunit = TUNIT_12HOURS;   break;
    case ISEC1_TABLE4_DAY:       timeunit = TUNIT_DAY;       break;
    default:
      if (lprint)
	{
	  Warning("GRIB time unit %d unsupported!", ISEC1_TimeUnit);
	  lprint = false;
	}
      break;
    }
  // clang-format on

  return timeunit;
}

static bool
cgribexTimeIsFC(int *isec1)
{
  bool isFC = (ISEC1_TimeRange == 10 && ISEC1_TimePeriod1 == 0 && ISEC1_TimePeriod2 == 0) ? false : true;
  return isFC;
}

static int
cgribexGetTsteptype(int timerange)
{
  static bool lprint = true;

  // clang-format off
  int tsteptype = TSTEP_INSTANT;
  switch ( timerange )
    {
    case  0:  tsteptype = TSTEP_INSTANT;  break;
    case  1:  tsteptype = TSTEP_INSTANT2; break;
    case  2:  tsteptype = TSTEP_RANGE;    break;
    case  3:  tsteptype = TSTEP_AVG;      break;
    case  4:  tsteptype = TSTEP_ACCUM;    break;
    case  5:  tsteptype = TSTEP_DIFF;     break;
    case 10:  tsteptype = TSTEP_INSTANT3; break;
    default:
      if (lprint)
	{
	  Warning("Time range indicator %d unsupported, set to 0!", timerange);
	  lprint = false;
	}
      break;
    }
  // clang-format on

  return tsteptype;
}

static bool
cgribexGetGridRegular(int *isec2, int *isec4, grid_t *grid, int gridtype, bool compyinc)
{
  bool ijDirectionIncrementGiven = gribbyte_get_bit(ISEC2_ResFlag, 1);
  bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);

  size_t nvalues = (size_t) ISEC4_NumValues;
  size_t nlon = (size_t) ISEC2_NumLon;
  size_t nlat = (size_t) ISEC2_NumLat;
  if (nvalues != nlon * nlat) Error("numberOfPoints (%zu) and gridSize (%zu) differ!", nvalues, nlon * nlat);

  grid->size = nvalues;
  grid->x.size = nlon;
  grid->y.size = nlat;

  if (gridtype == GRID_GAUSSIAN) grid->np = ISEC2_NumPar;
  grid->x.inc = 0;
  grid->y.inc = 0;
  grid->x.flag = 0;
  // if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 )
  {
    if (grid->x.size > 1)
      {
        bool recompinc = true;

        if (ISEC2_LastLon < ISEC2_FirstLon)
          {
            if (ISEC2_FirstLon >= 180000)
              ISEC2_FirstLon -= 360000;
            else
              ISEC2_LastLon += 360000;
          }

        if (ijDirectionIncrementGiven && ISEC2_LonIncr > 0)
          {
            if (labs(ISEC2_LastLon - (ISEC2_FirstLon + ISEC2_LonIncr * ((long) grid->x.size - 1))) <= 2)
              {
                recompinc = false;
                grid->x.inc = ISEC2_LonIncr * 0.001;
              }
          }

        // recompute xinc if necessary
        if (recompinc) grid->x.inc = (ISEC2_LastLon - ISEC2_FirstLon) * 0.001 / (grid->x.size - 1);

        // correct xinc if necessary
        if (ISEC2_FirstLon == 0 && ISEC2_LastLon > 354000 && ISEC2_LastLon < 360000)
          {
            double xinc = 360. / grid->x.size;
            if (fabs(grid->x.inc - xinc) > 0.0)
              {
                grid->x.inc = xinc;
                if (CDI_Debug) Message("set xinc to %g", grid->x.inc);
              }
          }
      }
    grid->x.first = ISEC2_FirstLon * 0.001;
    grid->x.last = ISEC2_LastLon * 0.001;
    grid->x.flag = 2;
  }
  grid->y.flag = 0;
  // if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 )
  {
    if (grid->y.size > 1 && compyinc)
      {
        bool recompinc = true;
        if (ijDirectionIncrementGiven && ISEC2_LatIncr > 0)
          {
            if (labs(ISEC2_LastLat - (ISEC2_FirstLat + ISEC2_LatIncr * ((long) grid->y.size - 1))) <= 2)
              {
                recompinc = false;
                grid->y.inc = ISEC2_LatIncr * 0.001;
              }
          }

        // recompute yinc if necessary
        if (recompinc) grid->y.inc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->y.size - 1);
      }
    grid->y.first = ISEC2_FirstLat * 0.001;
    grid->y.last = ISEC2_LastLat * 0.001;
    grid->y.flag = 2;
  }

  return uvRelativeToGrid;
}

static bool
cgribexGetGridReduced(int *isec2, int *isec4, grid_t *grid)
{
  bool ijDirectionIncrementGiven = gribbyte_get_bit(ISEC2_ResFlag, 1);
  bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);
  grid->np = ISEC2_NumPar;
  grid->size = (size_t) ISEC4_NumValues;

  size_t reducedPointsSize = (size_t) ISEC2_NumLat;
  grid->reducedPointsSize = reducedPointsSize;
  grid->reducedPoints = (int *) Malloc(reducedPointsSize * sizeof(int));
  memcpy(grid->reducedPoints, ISEC2_ReducedPointsPtr, reducedPointsSize * sizeof(int));

  grid->y.size = (size_t) ISEC2_NumLat;
  grid->x.inc = 0;
  grid->y.inc = 0;
  grid->x.flag = 0;
  // if ( ISEC2_FirstLon != 0 || ISEC2_LastLon != 0 )
  {
    if (ISEC2_LastLon < ISEC2_FirstLon)
      {
        if (ISEC2_FirstLon >= 180000)
          ISEC2_FirstLon -= 360000;
        else
          ISEC2_LastLon += 360000;
      }

    grid->x.first = ISEC2_FirstLon * 0.001;
    grid->x.last = ISEC2_LastLon * 0.001;
    grid->x.flag = 2;
  }
  grid->y.flag = 0;
  // if ( ISEC2_FirstLat != 0 || ISEC2_LastLat != 0 )
  {
    if (grid->y.size > 1)
      {
        if (ijDirectionIncrementGiven && ISEC2_LatIncr > 0)
          grid->y.inc = ISEC2_LatIncr * 0.001;
        else
          grid->y.inc = (ISEC2_LastLat - ISEC2_FirstLat) * 0.001 / (grid->y.size - 1);
      }
    grid->y.first = ISEC2_FirstLat * 0.001;
    grid->y.last = ISEC2_LastLat * 0.001;
    grid->y.flag = 2;
  }

  return uvRelativeToGrid;
}

static bool
cgribexGetGridLCC(int *isec2, int *isec4, grid_t *grid)
{
  bool uvRelativeToGrid = gribbyte_get_bit(ISEC2_ResFlag, 5);

  size_t nvalues = (size_t) ISEC4_NumValues;
  size_t nlon = (size_t) ISEC2_NumLon;
  size_t nlat = (size_t) ISEC2_NumLat;
  if (nvalues != nlon * nlat) Error("numberOfPoints (%zu) and gridSize (%zu) differ!", nvalues, nlon * nlat);

  grid->size = nvalues;
  grid->x.size = nlon;
  grid->y.size = nlat;

  grid->x.first = 0;
  grid->x.last = 0;
  grid->x.inc = ISEC2_Lambert_dx;
  grid->y.first = 0;
  grid->y.last = 0;
  grid->y.inc = ISEC2_Lambert_dy;
  grid->x.flag = 2;
  grid->y.flag = 2;

  return uvRelativeToGrid;
}

static bool
cgribexGetGrid(stream_t *streamptr, int *isec2, int *isec4, grid_t *grid, int iret)
{
  bool uvRelativeToGrid = false;
  bool compyinc = true;
  int gridtype = cgribexGetGridType(isec2);
  int projtype = (gridtype == GRID_PROJECTION && cgribexGetIsRotated(isec2)) ? CDI_PROJ_RLL : CDI_UNDEFID;
  if (gridtype == CDI_PROJ_LCC)
    {
      projtype = gridtype;
      gridtype = GRID_PROJECTION;
    }

  if (streamptr->unreduced && gridtype == GRID_GAUSSIAN_REDUCED && iret != -801)
    {
      int nlon = 0;
      for (int ilat = 0; ilat < ISEC2_NumLat; ++ilat)
        if (ISEC2_ReducedPoints(ilat) > nlon) nlon = ISEC2_ReducedPoints(ilat);
      gridtype = GRID_GAUSSIAN;
      ISEC2_NumLon = nlon;
      ISEC4_NumValues = nlon * ISEC2_NumLat;
      compyinc = false;
    }

  grid_init(grid);
  cdiGridTypeInit(grid, gridtype, 0);

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_RLL)
    {
      uvRelativeToGrid = cgribexGetGridRegular(isec2, isec4, grid, gridtype, compyinc);
    }
  else if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      uvRelativeToGrid = cgribexGetGridReduced(isec2, isec4, grid);
    }
  else if (projtype == CDI_PROJ_LCC)
    {
      uvRelativeToGrid = cgribexGetGridLCC(isec2, isec4, grid);
    }
  else if (gridtype == GRID_SPECTRAL)
    {
      grid->size = (size_t) ISEC4_NumValues;
      grid->trunc = ISEC2_PentaJ;
      grid->lcomplex = (ISEC2_RepMode == 2) ? 1 : 0;
    }
  else if (gridtype == GRID_GME)
    {
      grid->size = (size_t) ISEC4_NumValues;
      grid->gme.nd = ISEC2_GME_ND;
      grid->gme.ni = ISEC2_GME_NI;
      grid->gme.ni2 = ISEC2_GME_NI2;
      grid->gme.ni3 = ISEC2_GME_NI3;
    }
  else if (gridtype == GRID_GENERIC)
    {
      grid->size = (size_t) ISEC4_NumValues;
      grid->x.size = 0;
      grid->y.size = 0;
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }

  grid->type = gridtype;
  grid->projtype = projtype;

  return uvRelativeToGrid;
}

static void
cgribexGetLevel(int *isec1, int *leveltype, int *level1, int *level2)
{
  *leveltype = ISEC1_LevelType;
  *level1 = ISEC1_Level1;
  *level2 = ISEC1_Level2;
  if (*leveltype == GRIB1_LTYPE_ISOBARIC)
    *level1 *= 100;
  else if (*leveltype == GRIB1_LTYPE_99 || *leveltype == GRIB1_LTYPE_ISOBARIC_PA)
    *leveltype = GRIB1_LTYPE_ISOBARIC;
}

static void
cgribexDefProjLCC(int *isec2, int gridID)
{
  struct CDI_GridProjParams gpp;
  gridProjParamsInit(&gpp);

  bool earthIsOblate = gribbyte_get_bit(ISEC2_ResFlag, 2);
  if (earthIsOblate)
    {
      gpp.a = 6378160.0;
      gpp.b = 6356775.0;
      gpp.rf = 297.0;
    }
  else
    {
      gpp.a = 6367470.0;
    }

  gpp.xval_0 = ISEC2_FirstLon * 0.001;
  gpp.yval_0 = ISEC2_FirstLat * 0.001;
  gpp.lon_0 = ISEC2_Lambert_Lov * 0.001;
  gpp.lat_1 = ISEC2_Lambert_LatS1 * 0.001;
  gpp.lat_2 = ISEC2_Lambert_LatS2 * 0.001;
  bool lsouth = gribbyte_get_bit(ISEC2_Lambert_ProjFlag, 1);
  if (lsouth)
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
cgribexAddRecord(stream_t *streamptr, cgribexrec_t *cgribexp, int param, size_t recsize, off_t position, int comptype, int lmv,
                 int iret)
{
  int *isec1 = cgribexp->sec1;
  int *isec2 = cgribexp->sec2;
  int *isec4 = cgribexp->sec4;
  double *fsec2 = cgribexp->fsec2;
  double *fsec3 = cgribexp->fsec3;

  int datatype = (ISEC4_NumBits > 0 && ISEC4_NumBits <= 32) ? ISEC4_NumBits : CDI_DATATYPE_PACK;
  if (datatype > 32) datatype = CDI_DATATYPE_PACK32;
  if (datatype < 0) datatype = CDI_DATATYPE_PACK;

  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int recID = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  int tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);

  int leveltype, level1, level2;
  cgribexGetLevel(isec1, &leveltype, &level1, &level2);

  // fprintf(stderr, "param %d %d %d %d\n", param, level1, level2, leveltype);

  record->size = recsize;
  record->position = position;
  record->param = param;
  record->ilevel = level1;
  record->ilevel2 = level2;
  record->ltype = leveltype;
  record->tsteptype = (short) tsteptype;

  grid_t *gridptr = (grid_t *) Malloc(sizeof(*gridptr));
  bool uvRelativeToGrid = cgribexGetGrid(streamptr, isec2, isec4, gridptr, iret);

  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, gridptr, 0);
  int gridID = gridAdded.Id;
  if (!gridAdded.isNew)
    {
      grid_free(gridptr);
      Free(gridptr);
    }
  else if (gridptr->projtype == CDI_PROJ_RLL)
    {
      double xpole = ISEC2_LonSP * 0.001 - 180;
      double ypole = -ISEC2_LatSP * 0.001;
      double angle = -FSEC2_RotAngle;
      gridDefParamRLL(gridID, xpole, ypole, angle);
    }
  else if (gridptr->projtype == CDI_PROJ_LCC)
    {
      cgribexDefProjLCC(isec2, gridID);
    }

  int zaxistype = grib1ltypeToZaxisType(leveltype);
  if (zaxistype == ZAXIS_HYBRID || zaxistype == ZAXIS_HYBRID_HALF)
    {
      size_t vctsize = (size_t) ISEC2_NumVCP;
      double *vctptr = &fsec2[10];
      varDefVCT(vctsize, vctptr);
    }

  bool lbounds = cgribexGetZaxisHasBounds(leveltype);

  int varID = 0, levelID = 0;
  varAddRecord(recID, param, gridID, zaxistype, lbounds, level1, level2, 0, 0, datatype, &varID, &levelID, tsteptype, leveltype, -1,
               NULL, NULL, NULL, NULL);

  record->varID = (short) varID;
  record->levelID = levelID;

  varDefCompType(varID, comptype);

  if (uvRelativeToGrid) varDefKeyInt(varID, CDI_KEY_UVRELATIVETOGRID, 1);

  if (ISEC1_LocalFLag)
    {
      if (ISEC1_CenterID == 78 && isec1[36] == 253)  // DWD local extension
        {
          varDefKeyInt(varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, isec1[52]);
          varDefKeyInt(varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, isec1[53]);
          varDefKeyInt(varID, CDI_KEY_PERTURBATIONNUMBER, isec1[54]);
        }
      else if (ISEC1_CenterID == 252 && isec1[36] == 1)  // MPIM local extension
        {
          varDefKeyInt(varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, isec1[37]);
          varDefKeyInt(varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, isec1[39]);
          varDefKeyInt(varID, CDI_KEY_PERTURBATIONNUMBER, isec1[38]);
        }
    }

  if (lmv) varDefMissval(varID, FSEC3_MissVal);

  if (varInqInst(varID) == CDI_UNDEFID)
    {
      int center = ISEC1_CenterID;
      int subcenter = ISEC1_SubCenterID;
      int instID = institutInq(center, subcenter, NULL, NULL);
      if (instID == CDI_UNDEFID) instID = institutDef(center, subcenter, NULL, NULL);
      varDefInst(varID, instID);
    }

  if (varInqModel(varID) == CDI_UNDEFID)
    {
      int modelID = modelInq(varInqInst(varID), ISEC1_ModelID, NULL);
      if (modelID == CDI_UNDEFID) modelID = modelDef(varInqInst(varID), ISEC1_ModelID, NULL);
      varDefModel(varID, modelID);
    }

  if (varInqTable(varID) == CDI_UNDEFID)
    {
      int tableID = tableInq(varInqModel(varID), ISEC1_CodeTable, NULL);
      if (tableID == CDI_UNDEFID) tableID = tableDef(varInqModel(varID), ISEC1_CodeTable, NULL);
      varDefTable(varID, tableID);
    }

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
}

static void
MCH_get_undef(int *isec1, double *undef_pds, double *undef_eps)
{
  /* 2010-01-13: Oliver Fuhrer */
  if (ISEC1_CenterID == 215)
    {
      if (isec1[34] != 0 && isec1[34] != 255)
        {
          if (isec1[34] & 2)
            {
              *undef_pds = ((isec1[34] & 1) ? -0.99 : +0.99) * pow(10.0, -isec1[35]);
              *undef_eps = pow(10.0, -isec1[35] - 1);
            }
          else
            {
              *undef_pds = ((isec1[34] & 1) ? -0.99 : +0.99) * pow(10.0, +isec1[35]);
              *undef_eps = pow(10.0, isec1[35] - 1);
            }
        }
    }
}

static void
cgribexDecodeHeader(cgribexrec_t *cgribexp, int *gribbuffer, int recsize, int *lmv, int *iret)
{
  int *isec0 = cgribexp->sec0;
  int *isec1 = cgribexp->sec1;
  int *isec2 = cgribexp->sec2;
  int *isec3 = cgribexp->sec3;
  int *isec4 = cgribexp->sec4;
  double *fsec2 = cgribexp->fsec2;
  double *fsec3 = cgribexp->fsec3;

  int ipunp = 0, iword = 0;

  fill_intarr(isec1, 0, 256);
  fill_intarr(isec2, 0, 32);

  double *fsec4 = NULL;
  gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, fsec4, ipunp, (int *) gribbuffer, recsize, &iword, "J", iret);

  if (!(ISEC1_Sec2Or3Flag & 128)) isec2[0] = -1;  // default generic grid

  *lmv = 0;

  if (ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255))
    {
      double undef_pds, undef_eps;
      MCH_get_undef(isec1, &undef_pds, &undef_eps);
      FSEC3_MissVal = undef_pds;
      *lmv = 1;
    }
}

static compvar_t
cgribexVarSet(int param, int level1, int level2, int leveltype, int trange)
{
  int tsteptype = cgribexGetTsteptype(trange);

  compvar_t compVar;
  compVar.param = param;
  compVar.level1 = level1;
  compVar.level2 = level2;
  compVar.ltype = leveltype;
  compVar.tsteptype = tsteptype;

  return compVar;
}

static inline int
cgribexVarCompare(const compvar_t *compVar, const record_t *record, int flag)
{
  bool vinst
      = (compVar->tsteptype == TSTEP_INSTANT || compVar->tsteptype == TSTEP_INSTANT2 || compVar->tsteptype == TSTEP_INSTANT3);
  bool rinst = (record->tsteptype == TSTEP_INSTANT || record->tsteptype == TSTEP_INSTANT2 || record->tsteptype == TSTEP_INSTANT3);
  int tstepDiff = (!((flag == 0) & (vinst && rinst))) & (compVar->tsteptype != record->tsteptype);
  int rstatus = (compVar->param != record->param) | (compVar->level1 != record->ilevel) | (compVar->level2 != record->ilevel2)
                | (compVar->ltype != record->ltype) | tstepDiff;
  return rstatus;
}

#define gribWarning(text, nrecs, timestep, paramstr, level1, level2) \
  Warning("Record %2d (id=%s lev1=%d lev2=%d) timestep %d: %s", nrecs, paramstr, level1, level2, timestep, text)

static void
cgribexSkipRecords(int fileID)
{
  int nskip = CDI_Skip_Records;
  while (nskip-- > 0)
    {
      size_t recsize = gribGetSize(fileID);
      if (recsize == 0) Error("Skipping of %d records failed!", CDI_Skip_Records);

      off_t recpos = fileGetPos(fileID);
      fileSetPos(fileID, recpos, SEEK_CUR);
    }
}

static CdiDateTime
cgribexDateTimeX(int *isec1, CdiDateTime *sDateTime)
{
  int vdate = 0, sdate = 0, vtime = 0, stime = 0;
  gribDateTimeX(isec1, &vdate, &vtime, &sdate, &stime);

  sDateTime->date = cdiDate_set(sdate);
  sDateTime->time = cdiTime_set(stime);

  return cdiDateTime_set(vdate, vtime);
}

int
cgribexScanTimestep1(stream_t *streamptr)
{
  CdiDateTime vDateTime0;
  cdiDateTime_init(&vDateTime0);
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  void *gribbuffer = NULL;
  size_t buffersize = 0;
  int leveltype = 0, level1 = 0, level2 = 0;
  unsigned recID;
  int nrecsScanned = 0;
  bool warn_time = true;
  bool warn_numavg = true;
  bool fcast = false;
  char paramstr[32];

  streamptr->curTsID = 0;

  cgribexrec_t *cgribexp = (cgribexrec_t *) streamptr->record->objectp;

  int tsID = tstepsNewEntry(streamptr);
  if (tsID != 0) Error("Internal problem! tstepsNewEntry returns %d", tsID);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  int fileID = streamptr->fileID;

  if (CDI_Skip_Records) cgribexSkipRecords(fileID);

  unsigned nrecs = 0;
  while (true)
    {
      size_t recsize = gribGetSize(fileID);
      recpos = fileGetPos(fileID);

      if (recsize == 0)
        {
          if (nrecs == 0) Error("No GRIB records found!");
          streamptr->ntsteps = 1;
          break;
        }

      ensureBufferSize(recsize, &buffersize, &gribbuffer);

      size_t readsize = recsize;
      // Search for next 'GRIB', read the following record, and position file offset after it.
      if (gribRead(fileID, gribbuffer, &readsize)) break;

      int comptype = grbDecompress(recsize, &buffersize, &gribbuffer);

      size_t sec2len = cgribexSection2Length(gribbuffer, buffersize);
      if (sec2len > cgribexp->sec2len)
        {
          cgribexp->sec2len = sec2len;
          cgribexp->sec2 = (int *) Realloc(cgribexp->sec2, sec2len * sizeof(int));
        }

      int *isec1 = cgribexp->sec1;

      nrecsScanned++;
      cgribexDecodeHeader(cgribexp, (int *) gribbuffer, (int) recsize, &lmv, &iret);

      int param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      cgribexGetLevel(isec1, &leveltype, &level1, &level2);

      CdiDateTime sDateTime;
      CdiDateTime vDateTime = cgribexDateTimeX(isec1, &sDateTime);

      if (nrecs == 0)
        {
          vDateTime0 = vDateTime;
          fcast = cgribexTimeIsFC(isec1);
          taxis->unit = cgribexGetTimeUnit(isec1);
          taxis->rDateTime = cdiDateTime_set(gribRefDate(isec1), gribRefTime(isec1));
          taxis->sDateTime = sDateTime;
          taxis->vDateTime = vDateTime;
        }
      else
        {
          if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;

          compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);
          record_t *records = streamptr->tsteps[tsID].records;
          for (recID = 0; recID < nrecs; recID++)
            {
              if (cgribexVarCompare(&compVar, &records[recID], 0) == 0) break;
            }

          if (CDI_Inventory_Mode == 1)
            {
              if (recID < nrecs) break;
              if (warn_time)
                if (cdiDateTime_isNE(vDateTime, vDateTime0))
                  {
                    gribWarning("Inconsistent verification time!", nrecsScanned, tsID + 1, paramstr, level1, level2);
                    warn_time = false;
                  }
            }
          else
            {
              if (cdiDateTime_isNE(vDateTime, vDateTime0)) break;

              if (recID < nrecs)
                {
                  gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, paramstr, level1, level2);
                  continue;
                }
            }
        }

      if (ISEC1_AvgNum)
        {
          if (taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum))
            {
              Warning("Changing numavg from %d to %d not supported!", taxis->numavg, ISEC1_AvgNum);
              warn_numavg = false;
            }
          else
            {
              taxis->numavg = ISEC1_AvgNum;
            }
        }

      nrecs++;

      if (CDI_Debug)
        Message("Read record %2d (id=%s lev1=%d lev2=%d) %s", nrecsScanned, paramstr, level1, level2,
                CdiDateTime_string(vDateTime));

      cgribexAddRecord(streamptr, cgribexp, param, recsize, recpos, comptype, lmv, iret);
    }

  streamptr->rtsteps = 1;

  if (nrecs == 0) return CDI_EUFSTRUCT;

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
cgribexScanTimestep2(stream_t *streamptr)
{
  CdiDateTime vDateTime0;
  cdiDateTime_init(&vDateTime0);
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  int leveltype = 0, level1 = 0, level2 = 0;
  int recID = 0;
  bool warn_numavg = true;
  char paramstr[32];

  streamptr->curTsID = 1;

  cgribexrec_t *cgribexp = (cgribexrec_t *) streamptr->record->objectp;
  int *isec1 = cgribexp->sec1;
  int *isec2 = cgribexp->sec2;

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

  int nrecsScanned = nrecords;
  int rindex = 0;
  while (true)
    {
      if (rindex > nrecords) break;

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
      cgribexDecodeHeader(cgribexp, (int *) gribbuffer, (int) recsize, &lmv, &iret);

      int param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
      cdiParamToString(param, paramstr, sizeof(paramstr));

      cgribexGetLevel(isec1, &leveltype, &level1, &level2);

      CdiDateTime sDateTime;
      CdiDateTime vDateTime = cgribexDateTimeX(isec1, &sDateTime);

      if (rindex == 0)
        {
          vDateTime0 = vDateTime;
          int taxisID = vlistInqTaxis(vlistID);
          if (taxisInqType(taxisID) == TAXIS_RELATIVE)
            {
              taxis->type = TAXIS_RELATIVE;
              taxis->rDateTime = cdiDateTime_set(gribRefDate(isec1), gribRefTime(isec1));
            }
          else
            {
              taxis->type = TAXIS_ABSOLUTE;
            }
          taxis->unit = cgribexGetTimeUnit(isec1);
          taxis->vDateTime = vDateTime;
          taxis->sDateTime = sDateTime;
        }
      else
        {
          if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;
        }

      int tsteptype = cgribexGetTsteptype(ISEC1_TimeRange);

      if (ISEC1_AvgNum)
        {
          if (taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum))
            warn_numavg = false;
          else
            taxis->numavg = ISEC1_AvgNum;
        }

      compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);

      for (recID = 0; recID < nrecords; recID++)
        {
          if (cgribexVarCompare(&compVar, &records[recID], 0) == 0) break;
        }

      if (recID == nrecords)
        {
          gribWarning("Parameter not defined at timestep 1!", nrecsScanned, tsID + 1, paramstr, level1, level2);
          return CDI_EUFSTRUCT;
        }

      if (CDI_Inventory_Mode == 1)
        {
          if (records[recID].used)
            {
              break;
            }
          else
            {
              records[recID].used = true;
              streamptr->tsteps[tsID].recIDs[rindex] = recID;
            }
        }
      else
        {
          if (records[recID].used)
            {
              if (cdiDateTime_isNE(vDateTime, vDateTime0)) break;

              gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, paramstr, level1, level2);
              continue;
            }
          else
            {
              records[recID].used = true;
              streamptr->tsteps[tsID].recIDs[rindex] = recID;
            }
        }

      if (CDI_Debug)
        Message("Read record %2d (id=%s lev1=%d lev2=%d) %s", nrecsScanned, paramstr, level1, level2,
                CdiDateTime_string(vDateTime));

      if (cgribexVarCompare(&compVar, &streamptr->tsteps[tsID].records[recID], 0) != 0)
        {
          Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                  records[recID].ilevel, level1);
          return CDI_EUFSTRUCT;
        }

      records[recID].position = recpos;
      records[recID].size = recsize;

      int varID = records[recID].varID;
      int gridID = vlistInqVarGrid(vlistID, varID);
      if (gridInqSize(gridID) == 1 && gridInqType(gridID) == GRID_LONLAT)
        {
          if (IS_NOT_EQUAL(gridInqXval(gridID, 0), ISEC2_FirstLon * 0.001)
              || IS_NOT_EQUAL(gridInqYval(gridID, 0), ISEC2_FirstLat * 0.001))
            gridChangeType(gridID, GRID_TRAJECTORY);
        }

      if (tsteptype != TSTEP_INSTANT2 && tsteptype != vlistInqVarTsteptype(vlistID, varID))
        vlistDefVarTsteptype(vlistID, varID, tsteptype);

      rindex++;
    }

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

  return 0;
}

int
cgribexScanTimestep(stream_t *streamptr)
{
  CdiDateTime vDateTime0;
  cdiDateTime_init(&vDateTime0);
  int lmv = 0, iret = 0;
  off_t recpos = 0;
  int leveltype = 0, level1 = 0, level2 = 0;
  int vrecID, recID = 0;
  bool warn_numavg = true;
  int nrecs = 0;
  char paramstr[32];

  cgribexrec_t *cgribexp = (cgribexrec_t *) streamptr->record->objectp;
  int *isec1 = cgribexp->sec1;

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

      int nrecsScanned = streamptr->tsteps[0].nallrecs + streamptr->tsteps[1].nrecs * (tsID - 1);
      int rindex = 0;
      while (true)
        {
          if (rindex > nrecs) break;

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
          cgribexDecodeHeader(cgribexp, (int *) gribbuffer, (int) recsize, &lmv, &iret);

          int param = cdiEncodeParam(ISEC1_Parameter, ISEC1_CodeTable, 255);
          cdiParamToString(param, paramstr, sizeof(paramstr));

          cgribexGetLevel(isec1, &leveltype, &level1, &level2);

          CdiDateTime sDateTime;
          CdiDateTime vDateTime = cgribexDateTimeX(isec1, &sDateTime);

          if (rindex == nrecs) break;

          if (rindex == 0)
            {
              vDateTime0 = vDateTime;
              int vlistID = streamptr->vlistID;
              int taxisID = vlistInqTaxis(vlistID);
              if (taxisInqType(taxisID) == TAXIS_RELATIVE)
                {
                  taxis->type = TAXIS_RELATIVE;
                  taxis->rDateTime = cdiDateTime_set(gribRefDate(isec1), gribRefTime(isec1));
                }
              else
                {
                  taxis->type = TAXIS_ABSOLUTE;
                }
              taxis->unit = cgribexGetTimeUnit(isec1);
              taxis->vDateTime = vDateTime;
              taxis->sDateTime = sDateTime;
            }
          else
            {
              if (cdiDateTime_isLT(sDateTime, taxis->sDateTime)) taxis->sDateTime = sDateTime;
            }

          if (ISEC1_AvgNum)
            {
              if (taxis->numavg && warn_numavg && (taxis->numavg != ISEC1_AvgNum))
                warn_numavg = false;
              else
                taxis->numavg = ISEC1_AvgNum;
            }

          compvar_t compVar = cgribexVarSet(param, level1, level2, leveltype, ISEC1_TimeRange);

          for (vrecID = 0; vrecID < nrecs; vrecID++)
            {
              recID = streamptr->tsteps[1].recIDs[vrecID];
              if (cgribexVarCompare(&compVar, &records[recID], 0) == 0) break;
            }

          if (vrecID == nrecs)
            {
              gribWarning("Parameter not defined at timestep 1!", nrecsScanned, tsID + 1, paramstr, level1, level2);

              if (CDI_Inventory_Mode == 1)
                return CDI_EUFSTRUCT;
              else
                continue;
            }

          if (CDI_Inventory_Mode == 1)
            {
              records[recID].used = true;
              streamptr->tsteps[tsID].recIDs[rindex] = recID;
            }
          else
            {
              if (records[recID].used)
                {
                  char paramstr_[32];
                  cdiParamToString(param, paramstr_, sizeof(paramstr_));

                  if (cdiDateTime_isNE(vDateTime, vDateTime0)) break;

                  if (CDI_Debug)
                    gribWarning("Parameter already exist, skipped!", nrecsScanned, tsID + 1, paramstr_, level1, level2);

                  continue;
                }
              else
                {
                  records[recID].used = true;
                  streamptr->tsteps[tsID].recIDs[rindex] = recID;
                }
            }

          if (CDI_Debug)
            Message("Read record %2d (id=%s lev1=%d lev2=%d) %s", nrecsScanned, paramstr, level1, level2,
                    CdiDateTime_string(vDateTime));

          if (cgribexVarCompare(&compVar, &records[recID], 0) != 0)
            {
              Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                      records[recID].ilevel, level1);
              Error("Invalid, unsupported or inconsistent record structure");
            }

          records[recID].position = recpos;
          records[recID].size = recsize;

          rindex++;
        }

      for (vrecID = 0; vrecID < nrecs; vrecID++)
        {
          recID = streamptr->tsteps[tsID].recIDs[vrecID];
          if (!records[recID].used) break;
        }

      if (vrecID < nrecs)
        {
          cdiParamToString(records[recID].param, paramstr, sizeof(paramstr));
          gribWarning("Parameter not found!", nrecsScanned, tsID + 1, paramstr, records[recID].ilevel, records[recID].ilevel2);
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

int
cgribexDecode(int memtype, void *cgribex, void *gribbuffer, size_t gribsize, void *data, size_t datasize, int unreduced,
              size_t *nmiss, double missval)
{
  int status = 0;

  bool lalloc = cgribex == NULL;
  cgribexrec_t *cgribexp = (cgribexrec_t *) (lalloc ? cgribexNew() : cgribex);

  int *isec0 = cgribexp->sec0;
  int *isec1 = cgribexp->sec1;
  int *isec2 = cgribexp->sec2;
  int *isec3 = cgribexp->sec3;
  int *isec4 = cgribexp->sec4;
  double *fsec2 = cgribexp->fsec2;
  double *fsec3 = cgribexp->fsec3;
  float fsec2f[sizeof(cgribexp->fsec2) / sizeof(double)];
  float fsec3f[sizeof(cgribexp->fsec3) / sizeof(double)];

  char hoper[2];
  strcpy(hoper, unreduced ? "R" : "D");

  FSEC3_MissVal = missval;

  int iret = 0, iword = 0;
  if (memtype == MEMTYPE_FLOAT)
    gribExSP(isec0, isec1, isec2, fsec2f, isec3, fsec3f, isec4, (float *) data, (int) datasize, (int *) gribbuffer, (int) gribsize,
             &iword, hoper, &iret);
  else
    gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double *) data, (int) datasize, (int *) gribbuffer, (int) gribsize,
             &iword, hoper, &iret);

  *nmiss = (ISEC1_Sec2Or3Flag & 64) ? ISEC4_NumValues - ISEC4_NumNonMissValues : 0;

  if (ISEC1_CenterID == 215 && (isec1[34] != 0 && isec1[34] != 255))
    {
      double undef_pds, undef_eps;
      MCH_get_undef(isec1, &undef_pds, &undef_eps);

      *nmiss = 0;
      if (memtype == MEMTYPE_FLOAT)
        {
          float *restrict dataf = (float *) data;
          for (size_t i = 0; i < datasize; i++)
            if ((fabs(dataf[i] - undef_pds) < undef_eps) || IS_EQUAL(dataf[i], FSEC3_MissVal))
              {
                dataf[i] = (float) missval;
                (*nmiss)++;
              }
        }
      else
        {
          double *restrict datad = (double *) data;
          for (size_t i = 0; i < datasize; i++)
            if ((fabs(datad[i] - undef_pds) < undef_eps) || IS_EQUAL(datad[i], FSEC3_MissVal))
              {
                datad[i] = missval;
                (*nmiss)++;
              }
        }
    }

  if (lalloc) cgribexDelete(cgribexp);

  return status;
}

static void
cgribexDefInstitut(int *isec1, int vlistID, int varID)
{
  int instID = (vlistInqInstitut(vlistID) != CDI_UNDEFID) ? vlistInqInstitut(vlistID) : vlistInqVarInstitut(vlistID, varID);
  if (instID != CDI_UNDEFID)
    {
      ISEC1_CenterID = institutInqCenter(instID);
      ISEC1_SubCenterID = institutInqSubcenter(instID);
    }
}

static void
cgribexDefModel(int *isec1, int vlistID, int varID)
{
  int modelID = (vlistInqModel(vlistID) != CDI_UNDEFID) ? vlistInqModel(vlistID) : vlistInqVarModel(vlistID, varID);
  if (modelID != CDI_UNDEFID) ISEC1_ModelID = modelInqGribID(modelID);
}

static void
cgribexDefParam(int *isec1, int param)
{
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  if (pnum < 0) pnum = -pnum;

  static bool lwarn_pdis = true;
  if (pdis != 255 && lwarn_pdis)
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      Warning("Can't convert GRIB2 parameter ID (%s) to GRIB1, set to %d.%d!", paramstr, pnum, pcat);
      lwarn_pdis = false;
    }

  static bool lwarn_pnum = true;
  if (pnum > 255 && lwarn_pnum)
    {
      Warning("Parameter number %d out of range (1-255), set to %d!", pnum, pnum % 256);
      lwarn_pnum = false;
      pnum = pnum % 256;
    }

  ISEC1_CodeTable = pcat;
  ISEC1_Parameter = pnum;
}

static int
cgribexDefTimerange(int tsteptype, int factor, int calendar, CdiDateTime rDateTime, CdiDateTime vDateTime, CdiDateTime sDateTime,
                    int *pip1, int *pip2)
{
  JulianDate julianDate1 = julianDate_encode(calendar, rDateTime);
  JulianDate julianDate2 = julianDate_encode(calendar, vDateTime);
  JulianDate julianDate = julianDate_sub(julianDate2, julianDate1);

  int timerange = -1;
  int ip1 = 0, ip2 = 0;
  if (!(int) (fmod(julianDate_to_seconds(julianDate), factor)))
    {
      int ip = (int) lround(julianDate_to_seconds(julianDate) / factor);
      if ((ip > 255) && (tsteptype == TSTEP_INSTANT)) tsteptype = TSTEP_INSTANT3;

      int ipx = 0;
      if (!cdiDateTime_isNull(sDateTime)
          && (tsteptype == TSTEP_RANGE || tsteptype == TSTEP_AVG || tsteptype == TSTEP_ACCUM || tsteptype == TSTEP_DIFF))
        {
          julianDate2 = julianDate_encode(calendar, sDateTime);
          ipx = (int) lround(julianDate_to_seconds(julianDate_sub(julianDate2, julianDate1)) / factor);
        }

      // clang-format off
      switch (tsteptype)
	{
	case TSTEP_INSTANT:  timerange =  0; ip1 = ip;  ip2 = 0;  break;
	case TSTEP_INSTANT2: timerange =  1; ip1 = 0;   ip2 = 0;  break;
	case TSTEP_RANGE:    timerange =  2; ip1 = 0;   ip2 = ip; break;
	case TSTEP_AVG:      timerange =  3; ip1 = 0;   ip2 = ip; break;
	case TSTEP_ACCUM:    timerange =  4; ip1 = ipx; ip2 = ip; break;
	case TSTEP_DIFF:     timerange =  5; ip1 = 0;   ip2 = ip; break;
	case TSTEP_INSTANT3:
	default:             timerange = 10; ip1 = ip/256; ip2 = ip%256; break;
	}
      // clang-format on
    }

  *pip1 = ip1;
  *pip2 = ip2;

  return timerange;
}

static int
cgribexDefDateTime(int *isec1, int timeunit, CdiDateTime dt)
{
  int year, month, day, hour, minute, second, ms;
  cdiDate_decode(dt.date, &year, &month, &day);
  cdiTime_decode(dt.time, &hour, &minute, &second, &ms);

  int century = year / 100;
  ISEC1_Year = year - century * 100;

  if (year < 0)
    {
      century = -century;
      ISEC1_Year = -ISEC1_Year;
    }

  if (ISEC1_Year == 0)
    {
      century -= 1;
      ISEC1_Year = 100;
    }

  century += 1;
  if (year < 0) century = -century;

  ISEC1_Month = month;
  ISEC1_Day = day;
  ISEC1_Hour = hour;
  ISEC1_Minute = minute;

  ISEC1_Century = century;

  int factor = 1;
  // clang-format off
  switch (timeunit)
    {
    case TUNIT_MINUTE:    factor =    60; ISEC1_TimeUnit = ISEC1_TABLE4_MINUTE;    break;
    case TUNIT_QUARTER:   factor =   900; ISEC1_TimeUnit = ISEC1_TABLE4_QUARTER;   break;
    case TUNIT_30MINUTES: factor =  1800; ISEC1_TimeUnit = ISEC1_TABLE4_30MINUTES; break;
    case TUNIT_HOUR:      factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    case TUNIT_3HOURS:    factor = 10800; ISEC1_TimeUnit = ISEC1_TABLE4_3HOURS;    break;
    case TUNIT_6HOURS:    factor = 21600; ISEC1_TimeUnit = ISEC1_TABLE4_6HOURS;    break;
    case TUNIT_12HOURS:   factor = 43200; ISEC1_TimeUnit = ISEC1_TABLE4_12HOURS;   break;
    case TUNIT_DAY:       factor = 86400; ISEC1_TimeUnit = ISEC1_TABLE4_DAY;       break;
    default:              factor =  3600; ISEC1_TimeUnit = ISEC1_TABLE4_HOUR;      break;
    }
  // clang-format on

  return factor;
}

static void
cgribexDefTime(int *isec1, CdiDateTime vDateTime, int tsteptype, int numavg, int taxisID)
{
  int timetype = TAXIS_ABSOLUTE;
  int timeunit = TUNIT_HOUR;

  if (taxisID != -1)
    {
      timetype = taxisInqType(taxisID);
      timeunit = taxisInqTunit(taxisID);
    }

  if (timetype == TAXIS_RELATIVE)
    {
      int ip1 = 0, ip2 = 0;
      int calendar = taxisInqCalendar(taxisID);

      CdiDateTime rDateTime = taxisInqRdatetime(taxisID);
      if (cdiDateTime_isLT(vDateTime, rDateTime)) rDateTime = vDateTime;

      CdiDateTime sDateTime = taxisInqSdatetime(taxisID);

      int factor = cgribexDefDateTime(isec1, timeunit, rDateTime);
      int timerange = cgribexDefTimerange(tsteptype, factor, calendar, rDateTime, vDateTime, sDateTime, &ip1, &ip2);

      if (ip2 > 0xFF)
        {
          rDateTime = vDateTime;
          factor = cgribexDefDateTime(isec1, timeunit, rDateTime);
          timerange = cgribexDefTimerange(tsteptype, factor, calendar, rDateTime, vDateTime, sDateTime, &ip1, &ip2);
        }
      /*
      if (ip2 > 0xFF && timeunit < TUNIT_YEAR)
        {
          timeunit++;
          factor = cgribexDefDateTime(isec1, timeunit, rDateTime);
          timerange = cgribexDefTimerange(tsteptype, factor, calendar, rDateTime, vDateTime, sDateTime, &ip1, &ip2);
        }
      */
      if (timerange == -1 || timerange == 1 || timerange == 3) timetype = TAXIS_ABSOLUTE;
      /*
      else if (timerange == 10)
        {
          if (ip1 < 0 || ip1 > 0xFFFF) timetype = TAXIS_ABSOLUTE;
          if (ip2 < 0 || ip2 > 0xFFFF) timetype = TAXIS_ABSOLUTE;
        }
      */
      else
        {
          if (ip1 < 0 || ip1 > 0xFF) timetype = TAXIS_ABSOLUTE;
          if (ip2 < 0 || ip2 > 0xFF) timetype = TAXIS_ABSOLUTE;
        }

      if (timetype != TAXIS_ABSOLUTE)
        {
          ISEC1_TimeRange = timerange;
          ISEC1_TimePeriod1 = ip1;
          ISEC1_TimePeriod2 = ip2;
        }
    }

  if (timetype == TAXIS_ABSOLUTE)
    {
      (void) cgribexDefDateTime(isec1, timeunit, vDateTime);

      /*
      if (numavg > 0)
        ISEC1_TimeRange = 0;
      else
      */
      if (ISEC1_TimeRange != 3) ISEC1_TimeRange = 10;

      ISEC1_TimePeriod1 = 0;
      ISEC1_TimePeriod2 = 0;
    }

  ISEC1_AvgNum = numavg;
  ISEC1_AvgMiss = 0;
  ISEC1_DecScaleFactor = 0;
}

static void
cgribexDefGridRegular(int *isec2, double *fsec2, int gridID, int gridtype, bool gridIsRotated, bool gridIsCurvilinear,
                      int uvRelativeToGrid)
{
  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    ISEC2_GridType = GRIB1_GTYPE_GAUSSIAN;
  else if (gridtype == GRID_LONLAT && gridIsRotated)
    ISEC2_GridType = GRIB1_GTYPE_LATLON_ROT;
  else
    ISEC2_GridType = GRIB1_GTYPE_LATLON;

  double xfirst = 0.0, xlast = 0.0, xinc = 0.0;
  double yfirst = 0.0, ylast = 0.0, yinc = 0.0;

  int nlon = (int) gridInqXsize(gridID);
  int nlat = (int) gridInqYsize(gridID);

  if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      ISEC2_Reduced = true;
      if (nlon == 2)
        {
          xfirst = gridInqXval(gridID, 0);
          xlast = gridInqXval(gridID, 1);
        }
      else
        {
          xlast = 360.0 - 360.0 / (nlat * 2);
        }

      nlon = 0;
      gridInqReducedPoints(gridID, ISEC2_ReducedPointsPtr);
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

  ISEC2_NumLon = nlon;
  ISEC2_NumLat = nlat;
  ISEC2_FirstLat = (int) lround(yfirst * 1000);
  ISEC2_LastLat = (int) lround(ylast * 1000);
  ISEC2_FirstLon = (int) lround(xfirst * 1000);
  ISEC2_LastLon = (int) lround(xlast * 1000);
  // gribapi gridType detector doesn't like lonIncr for Gaussian reduced longitides
  if (gridtype != GRID_GAUSSIAN_REDUCED) ISEC2_LonIncr = (int) lround(xinc * 1000);

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      int np = gridInqNP(gridID);
      if (np == 0) np = nlat / 2;
      ISEC2_NumPar = np;
    }
  else
    {
      ISEC2_LatIncr = (int) lround(yinc * 1000);
    }

  if (ISEC2_NumLon > 1 && ISEC2_NumLat == 1)
    if (ISEC2_LonIncr != 0 && ISEC2_LatIncr == 0) ISEC2_LatIncr = ISEC2_LonIncr;

  if (ISEC2_NumLon == 1 && ISEC2_NumLat > 1)
    if (ISEC2_LonIncr == 0 && ISEC2_LatIncr != 0) ISEC2_LonIncr = ISEC2_LatIncr;

  ISEC2_ResFlag = 0;
  if (ISEC2_LatIncr && ISEC2_LonIncr) gribbyte_set_bit(&ISEC2_ResFlag, 1);
  if (uvRelativeToGrid > 0) gribbyte_set_bit(&ISEC2_ResFlag, 5);

  if (gridIsRotated)
    {
      double xpole = 0, ypole = 0, angle = 0;
      gridInqParamRLL(gridID, &xpole, &ypole, &angle);

      ISEC2_LatSP = -(int) lround(ypole * 1000);
      ISEC2_LonSP = (int) lround((xpole + 180) * 1000);
      if (fabs(angle) > 0) angle = -angle;
      FSEC2_RotAngle = angle;
    }

  ISEC2_ScanFlag = 0;
  if (ISEC2_LastLon < ISEC2_FirstLon) gribbyte_set_bit(&ISEC2_ScanFlag, 1);  // East -> West
  if (ISEC2_LastLat > ISEC2_FirstLat) gribbyte_set_bit(&ISEC2_ScanFlag, 2);  // South -> North
}

static void
cgribexDefGridLambert(int *isec2, int gridID, int uvRelativeToGrid)
{
  int xsize = (int) gridInqXsize(gridID);
  int ysize = (int) gridInqYsize(gridID);

  struct CDI_GridProjParams gpp;
  gridInqParamsLCC(gridID, &gpp);
  if (IS_EQUAL(gpp.x_0, gpp.mv) && IS_EQUAL(gpp.y_0, gpp.mv) && (IS_EQUAL(gpp.xval_0, gpp.mv) || IS_EQUAL(gpp.yval_0, gpp.mv)))
    {
      gpp.x_0 = gridInqXval(gridID, 0);
      gpp.y_0 = gridInqYval(gridID, 0);
    }
  gridVerifyProjParamsLCC(&gpp);

  bool lsouth = (gpp.lat_1 < 0);
  if (lsouth)
    {
      gpp.lat_1 = -gpp.lat_2;
      gpp.lat_2 = -gpp.lat_2;
    }

  double xinc = gridInqXinc(gridID);
  double yinc = gridInqYinc(gridID);
  if (IS_EQUAL(xinc, 0.0)) xinc = gridInqXincInMeter(gridID);
  if (IS_EQUAL(yinc, 0.0)) yinc = gridInqYincInMeter(gridID);

  ISEC2_GridType = GRIB1_GTYPE_LCC;
  ISEC2_NumLon = xsize;
  ISEC2_NumLat = ysize;
  ISEC2_FirstLon = (int) lround(gpp.xval_0 * 1000);
  ISEC2_FirstLat = (int) lround(gpp.yval_0 * 1000);
  ISEC2_Lambert_Lov = (int) lround(gpp.lon_0 * 1000);
  ISEC2_Lambert_LatS1 = (int) lround(gpp.lat_1 * 1000);
  ISEC2_Lambert_LatS2 = (int) lround(gpp.lat_2 * 1000);
  ISEC2_Lambert_dx = (int) lround(xinc);
  ISEC2_Lambert_dy = (int) lround(yinc);
  ISEC2_Lambert_LatSP = 0;
  ISEC2_Lambert_LonSP = 0;
  ISEC2_Lambert_ProjFlag = 0;
  if (lsouth) gribbyte_set_bit(&ISEC2_Lambert_ProjFlag, 1);

  bool earthIsOblate = (IS_EQUAL(gpp.a, 6378160.0) && IS_EQUAL(gpp.rf, 297.0));
  ISEC2_ResFlag = 0;
  if (ISEC2_Lambert_dx && ISEC2_Lambert_dy) gribbyte_set_bit(&ISEC2_ResFlag, 1);
  if (earthIsOblate) gribbyte_set_bit(&ISEC2_ResFlag, 2);
  if (uvRelativeToGrid > 0) gribbyte_set_bit(&ISEC2_ResFlag, 5);

  ISEC2_ScanFlag = 0;
  gribbyte_set_bit(&ISEC2_ScanFlag, 2);  // South -> North
}

static void
cgribexDefGridSpectal(int *isec2, int *isec4, int gridID)
{
  ISEC2_GridType = GRIB1_GTYPE_SPECTRAL;
  ISEC2_PentaJ = gridInqTrunc(gridID);
  ISEC2_PentaK = ISEC2_PentaJ;
  ISEC2_PentaM = ISEC2_PentaJ;
  ISEC2_RepType = 1;
  isec4[2] = 128;
  if (gridInqComplexPacking(gridID) && ISEC2_PentaJ >= 21)
    {
      ISEC2_RepMode = 2;
      isec4[3] = 64;
      isec4[16] = 0;
      isec4[17] = 20;
      isec4[18] = 20;
      isec4[19] = 20;
    }
  else
    {
      ISEC2_RepMode = 1;
      isec4[3] = 0;
    }
}

static void
cgribexDefGridGME(int *isec2, int gridID)
{
  ISEC2_GridType = GRIB1_GTYPE_GME;
  int nd = 0, ni = 0, ni2 = 0, ni3 = 0;
  gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
  ISEC2_GME_ND = nd;
  ISEC2_GME_NI = ni;
  ISEC2_GME_NI2 = ni2;
  ISEC2_GME_NI3 = ni3;
  ISEC2_GME_AFlag = 0;
  ISEC2_GME_LatPP = 90000;
  ISEC2_GME_LonPP = 0;
  ISEC2_GME_LonMPL = 0;
  ISEC2_GME_BFlag = 0;
}

static void
cgribexDefGrid(int *isec1, int *isec2, double *fsec2, int *isec4, int gridID, int uvRelativeToGrid)
{
  fill_intarr(isec2, 0, 16);
  ISEC1_Sec2Or3Flag = 128;
  ISEC1_GridDefinition = 255;
  ISEC2_Reduced = false;
  ISEC2_ScanFlag = 0;

  int gridsize = (int) gridInqSize(gridID);
  bool gridIsRotated = false;
  bool gridIsCurvilinear = false;
  int gridtype = grbGetGridtype(&gridID, gridsize, &gridIsRotated, &gridIsCurvilinear);

  switch (gridtype)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_TRAJECTORY:
      {
        cgribexDefGridRegular(isec2, fsec2, gridID, gridtype, gridIsRotated, gridIsCurvilinear, uvRelativeToGrid);
        break;
      }
    case CDI_PROJ_LCC:
      {
        cgribexDefGridLambert(isec2, gridID, uvRelativeToGrid);
        break;
      }
    case GRID_SPECTRAL:
      {
        cgribexDefGridSpectal(isec2, isec4, gridID);
        break;
      }
    case GRID_GME:
      {
        cgribexDefGridGME(isec2, gridID);
        break;
      }
    case GRID_GENERIC:
      {
        ISEC1_Sec2Or3Flag = 0;
        break;
      }
    default:
      {
        static bool lwarn = true;
        ISEC1_Sec2Or3Flag = 0;
        if (lwarn) Warning("CGRIBEX library doesn't support %s grids, grid information will be lost!", gridNamePtr(gridtype));
        lwarn = false;
        break;
      }
    }
}

static int
level2int(double level)
{
  return (int) round(level);
}

static void
isec1DefLevel(int *isec1, int leveltype, int level1, int level2)
{
  ISEC1_LevelType = leveltype;
  ISEC1_Level1 = level1;
  ISEC1_Level2 = level2;
}

static void
cgribexDefLevel(int *isec1, int *isec2, double *fsec2, int zaxisID, int levelID)
{
  int zaxistype = zaxisInqType(zaxisID);
  int ltype = 0;
  cdiInqKeyInt(zaxisID, CDI_GLOBAL, CDI_KEY_TYPEOFFIRSTFIXEDSURFACE, &ltype);

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

  ISEC2_NumVCP = 0;

  int grib_ltype = zaxisTypeToGrib1ltype(zaxistype);

  switch (zaxistype)
    {
    case ZAXIS_SURFACE:
    case ZAXIS_MEANSEA:
    case ZAXIS_ALTITUDE:
    case ZAXIS_DEPTH_BELOW_SEA:
    case ZAXIS_ISENTROPIC:
      {
        isec1DefLevel(isec1, grib_ltype, level2int(level), 0);
        break;
      }
    case ZAXIS_CLOUD_BASE:
    case ZAXIS_CLOUD_TOP:
    case ZAXIS_ISOTHERM_ZERO:
    case ZAXIS_TROPOPAUSE:
    case ZAXIS_TOA:
    case ZAXIS_SEA_BOTTOM:
    case ZAXIS_ATMOSPHERE:
      {
        isec1DefLevel(isec1, grib_ltype, 0, 0);
        break;
      }
    case ZAXIS_HYBRID:
    case ZAXIS_HYBRID_HALF:
      {
        grib_ltype = hasBounds ? GRIB1_LTYPE_HYBRID_LAYER : GRIB1_LTYPE_HYBRID;
        isec1DefLevel(isec1, grib_ltype, level2int(dlevel1), level2int(dlevel2));

        int vctsize = zaxisInqVctSize(zaxisID);
        if (vctsize > 255)
          {
            static bool lwarning_vct = true;
            ISEC2_NumVCP = 0;
            if (lwarning_vct)
              {
                Warning("VCT size of %d is too large (maximum is 255). Set to 0!", vctsize);
                lwarning_vct = false;
              }
          }
        else
          {
            ISEC2_NumVCP = vctsize;
            zaxisInqVct(zaxisID, &fsec2[10]);
          }
        break;
      }
    case ZAXIS_PRESSURE:
      {
        if (level < 0) Warning("Pressure level of %f Pa is below zero!", level);

        if (!zaxis_units_is_Pa(zaxisID)) level *= 100.0;

        double dum;
        if (level < 32768 && (level < 100 || modf(level / 100, &dum) > 0))
          grib_ltype = GRIB1_LTYPE_ISOBARIC_PA;
        else
          level = level / 100;

        isec1DefLevel(isec1, grib_ltype, level2int(level), 0);
        break;
      }
    case ZAXIS_HEIGHT:
      {
        double sf = zaxis_units_to_meter(zaxisID);
        isec1DefLevel(isec1, grib_ltype, level2int(level * sf), 0);
        break;
      }
    case ZAXIS_SIGMA:
      {
        grib_ltype = hasBounds ? GRIB1_LTYPE_SIGMA_LAYER : GRIB1_LTYPE_SIGMA;
        isec1DefLevel(isec1, grib_ltype, level2int(dlevel1), level2int(dlevel2));
        break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
        grib_ltype = hasBounds ? GRIB1_LTYPE_LANDDEPTH_LAYER : GRIB1_LTYPE_LANDDEPTH;
        double sf = zaxis_units_to_centimeter(zaxisID);
        isec1DefLevel(isec1, grib_ltype, level2int(sf * dlevel1), level2int(sf * dlevel2));
        break;
      }
    case ZAXIS_GENERIC:
      {
        isec1DefLevel(isec1, ltype, level2int(level), 0);
        break;
      }
    default:
      {
        Error("Unsupported zaxis type: %s", zaxisNamePtr(zaxistype));
        break;
      }
    }
}

static void
cgribexDefaultSec0(int *isec0)
{
  ISEC0_GRIB_Len = 0;
  ISEC0_GRIB_Version = 0;
}

static void
cgribexDefaultSec1(int *isec1)
{
  ISEC1_CenterID = 0;
  ISEC1_SubCenterID = 0;
  ISEC1_LocalFLag = 0;
}

static void
cgribexDefaultSec4(int *isec4)
{
  for (int i = 2; i <= 10; ++i) isec4[i] = 0;
}

static void
cgribexDefEnsembleVar(int *isec1, int vlistID, int varID)
{
  // For Ensemble info

  // Put1Byte(isec1[36]);        // MPIM local GRIB use definition identifier (extension identifier)
  // Put1Byte(isec1[37]);        // type of ensemble forecast
  // Put2Byte(isec1[38]);        // individual ensemble member
  // Put2Byte(isec1[39]);        // number of forecasts in ensemble

  int perturbationNumber, numberOfForecastsInEnsemble, typeOfEnsembleForecast;
  int r1 = cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber);
  int r2 = cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble);
  int r3 = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, &typeOfEnsembleForecast);

  if (r1 == 0 && r2 == 0 && r3 == 0)
    {
      if (ISEC1_CenterID == 252)
        {
          ISEC1_LocalFLag = 1;
          isec1[36] = 1;
          isec1[37] = typeOfEnsembleForecast;
          isec1[38] = perturbationNumber;
          isec1[39] = numberOfForecastsInEnsemble;
        }
    }
}

size_t
cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID, CdiDateTime vDateTime, int tsteptype,
              int numavg, size_t datasize, const void *data, size_t nmiss, void *gribbuffer, size_t gribbuffersize)
{
  cgribexrec_t *cgribexp = (cgribexrec_t *) cgribexNew();

  size_t sec2len = 1024 + 2 * gridInqYsize(gridID);  // Gaussian reduced grid
  if (sec2len > cgribexp->sec2len)
    {
      cgribexp->sec2len = sec2len;
      cgribexp->sec2 = (int *) Realloc(cgribexp->sec2, sec2len * sizeof(int));
    }

  int *isec0 = cgribexp->sec0;
  int *isec1 = cgribexp->sec1;
  int *isec2 = cgribexp->sec2;
  int *isec3 = cgribexp->sec3;
  int *isec4 = cgribexp->sec4;
  double *fsec2 = cgribexp->fsec2;
  double *fsec3 = cgribexp->fsec3;
  float fsec2f[sizeof(cgribexp->fsec2) / sizeof(double)];
  float fsec3f[sizeof(cgribexp->fsec3) / sizeof(double)];

  fill_intarr(isec1, 0, 256);
  fsec2[0] = 0;
  fsec2[1] = 0;
  fsec2f[0] = 0;
  fsec2f[1] = 0;

  int gribsize = (int) (gribbuffersize / sizeof(int));
  int param = vlistInqVarParam(vlistID, varID);

  cgribexDefaultSec0(isec0);
  cgribexDefaultSec1(isec1);
  cgribexDefaultSec4(isec4);

  cgribexDefInstitut(isec1, vlistID, varID);
  cgribexDefModel(isec1, vlistID, varID);

  int datatype = vlistInqVarDatatype(vlistID, varID);

  int uvRelativeToGrid = -1;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_UVRELATIVETOGRID, &uvRelativeToGrid);

  cgribexDefParam(isec1, param);
  cgribexDefTime(isec1, vDateTime, tsteptype, numavg, vlistInqTaxis(vlistID));
  cgribexDefGrid(isec1, isec2, fsec2, isec4, gridID, uvRelativeToGrid);
  cgribexDefLevel(isec1, isec2, fsec2, zaxisID, levelID);

  cgribexDefEnsembleVar(isec1, vlistID, varID);

  cdi_check_gridsize_int_limit("GRIB1", datasize);

  ISEC4_NumValues = (int) datasize;
  ISEC4_NumBits = grbBitsPerValue(datatype);

  if (nmiss > 0)
    {
      FSEC3_MissVal = vlistInqVarMissval(vlistID, varID);
      ISEC1_Sec2Or3Flag |= 64;
    }

  if (isec4[2] == 128 && isec4[3] == 64)
    {
      if (memtype == MEMTYPE_FLOAT)
        isec4[16] = (int) (1000 * calculate_pfactor_float((const float *) data, ISEC2_PentaJ, isec4[17]));
      else
        isec4[16] = (int) (1000 * calculate_pfactor_double((const double *) data, ISEC2_PentaJ, isec4[17]));
      if (isec4[16] < -10000) isec4[16] = -10000;
      if (isec4[16] > 10000) isec4[16] = 10000;
    }
  // printf("isec4[16] %d\n", isec4[16]);

  if (memtype == MEMTYPE_FLOAT)
    {
      int numVCP = (ISEC2_NumVCP > 0) ? ISEC2_NumVCP : 0;
      for (int i = 0; i < numVCP; ++i) fsec2f[10 + i] = (float) fsec2[10 + i];
      fsec3f[1] = (float) fsec3[1];
    }

  int iret = 0, iword = 0;
  if (memtype == MEMTYPE_FLOAT)
    gribExSP(isec0, isec1, isec2, fsec2f, isec3, fsec3f, isec4, (float *) data, (int) datasize, (int *) gribbuffer, gribsize,
             &iword, "C", &iret);
  else
    gribExDP(isec0, isec1, isec2, fsec2, isec3, fsec3, isec4, (double *) data, (int) datasize, (int *) gribbuffer, gribsize, &iword,
             "C", &iret);

  cgribexDelete(cgribexp);

  if (iret) Error("Problem during GRIB encode (errno = %d)!", iret);

  size_t nbytes = (size_t) iword * sizeof(int);
  return nbytes;
}

void
cgribexChangeParameterIdentification(void *gh, int code, int ltype, int lev)
{
  if (!gh) return;

  unsigned char *pds = ((cgribex_handle *) gh)->pds;
  if (!pds) return;

  pds[8] = (unsigned char) code;
  pds[9] = (unsigned char) ltype;
  pds[10] = (unsigned char) lev;
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
