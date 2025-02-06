#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <float.h>

#include "dmemory.h"

#include "error.h"
#include "file.h"
#include "cdi.h"
#include "cdi_int.h"
#include "varscan.h"
#include "ieg.h"
#include "stream_scan.h"
#include "stream_ieg.h"
#include "get_num_missvals.h"
#include "exse.h"

#ifdef HAVE_LIBIEG

static int
iegInqDatatype(int prec)
{
  return (prec == EXSE_DOUBLE_PRECISION) ? CDI_DATATYPE_FLT64 : CDI_DATATYPE_FLT32;
}

static int
iegDefDatatype(int datatype)
{
  if (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64) Error("CDI/IEG library does not support complex numbers!");

  if (datatype != CDI_DATATYPE_FLT32 && datatype != CDI_DATATYPE_FLT64) datatype = CDI_DATATYPE_FLT32;

  return (datatype == CDI_DATATYPE_FLT64) ? EXSE_DOUBLE_PRECISION : EXSE_SINGLE_PRECISION;
}

static void
ieg_read_recordSP(stream_t *streamptr, float *data, size_t *numMissVals)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;

  int vrecID = streamptr->tsteps[tsID].curRecID;
  int recID = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *iegp = streamptr->record->objectp;
  if (iegRead(fileID, iegp) < 0) Error("Could not read IEG record!");

  iegInqDataSP(iegp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t size = gridInqSize(vlistInqVarGrid(vlistID, varID));

  *numMissVals = get_num_missvalsSP(size, data, missval);

  streamptr->numvals += size;
}

static void
ieg_read_recordDP(stream_t *streamptr, double *data, size_t *numMissVals)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;

  int vrecID = streamptr->tsteps[tsID].curRecID;
  int recID = streamptr->tsteps[tsID].recIDs[vrecID];
  int varID = streamptr->tsteps[tsID].records[recID].varID;
  off_t recpos = streamptr->tsteps[tsID].records[recID].position;

  fileSetPos(fileID, recpos, SEEK_SET);

  void *iegp = streamptr->record->objectp;
  if (iegRead(fileID, iegp) < 0) Error("Could not read IEG record!");

  iegInqDataDP(iegp, data);

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t size = gridInqSize(vlistInqVarGrid(vlistID, varID));

  *numMissVals = get_num_missvalsDP(size, data, missval);

  streamptr->numvals += size;
}

void
ieg_read_record(stream_t *streamptr, int memtype, void *data, size_t *numMissVals)
{
  if (memtype == MEMTYPE_DOUBLE)
    ieg_read_recordDP(streamptr, (double *) data, numMissVals);
  else
    ieg_read_recordSP(streamptr, (float *) data, numMissVals);
}

static int
iegGetZaxisType(int iegleveltype)
{
  // clang-format off
  switch (iegleveltype)
    {
    case IEG_LTYPE_SURFACE:         return ZAXIS_SURFACE;
    case IEG_LTYPE_99:
    case IEG_LTYPE_ISOBARIC:        return ZAXIS_PRESSURE;
    case IEG_LTYPE_HEIGHT:          return ZAXIS_HEIGHT;
    case IEG_LTYPE_ALTITUDE:        return ZAXIS_ALTITUDE;
    case IEG_LTYPE_HYBRID:
    case IEG_LTYPE_HYBRID_LAYER:    return ZAXIS_HYBRID;
    case IEG_LTYPE_LANDDEPTH:
    case IEG_LTYPE_LANDDEPTH_LAYER: return ZAXIS_DEPTH_BELOW_LAND;
    case IEG_LTYPE_SEADEPTH:        return ZAXIS_DEPTH_BELOW_SEA;
    }
  // clang-format on
  return ZAXIS_GENERIC;
}

static void
iegDefTime(int *pdb, CdiDateTime datetime, int taxisID)
{
  int timetype = (taxisID != -1) ? taxisInqType(taxisID) : -1;
  if (timetype == TAXIS_ABSOLUTE || timetype == TAXIS_RELATIVE)
    {
      IEG_P_Year(pdb) = datetime.date.year;
      IEG_P_Month(pdb) = datetime.date.month;
      IEG_P_Day(pdb) = datetime.date.day;
      IEG_P_Hour(pdb) = datetime.time.hour;
      IEG_P_Minute(pdb) = datetime.time.minute;

      pdb[15] = 1;
      pdb[16] = 0;
      pdb[17] = 0;
      pdb[18] = 10;
      pdb[36] = 1;
    }

  pdb[5] = 128;
}

/* find smallest power of 10 in [1000,10000000] that upon multiplication results in fractional part
   close to zero for all arguments */
static double
calc_resfac(double xfirst, double xlast, double xinc, double yfirst, double ylast, double yinc)
{
  double resfac = 1000.0;
  enum
  {
    nPwrOf10 = 5,
    nMultTests = 6,
  };
  static const double scaleFactors[nPwrOf10] = { 1000, 10000, 100000, 1000000, 10000000 };
  const double vals[nMultTests] = { xfirst, xlast, xinc, yfirst, ylast, yinc };

  for (size_t j = 0; j < nPwrOf10; ++j)
    {
      double scaleBy = scaleFactors[j];
      bool fractionalScale = false;
      for (size_t i = 0; i < nMultTests; ++i)
        {
          fractionalScale = fractionalScale || fabs(vals[i] * scaleBy - round(vals[i] * scaleBy)) > FLT_EPSILON;
        }
      if (!fractionalScale)
        {
          resfac = scaleBy;
          break;
        }
    }

  return resfac;
}

static void
iegDefGrid(int *gdb, int gridID)
{
  int projID = gridInqProj(gridID);
  if (projID != CDI_UNDEFID && gridInqProjType(projID) == CDI_PROJ_RLL) gridID = projID;

  int gridtype = gridInqType(gridID);
  int projtype = CDI_UNDEFID;
  if (gridtype == GRID_PROJECTION && gridInqProjType(gridID) == CDI_PROJ_RLL) projtype = CDI_PROJ_RLL;

  size_t xsize = gridInqXsize(gridID);
  size_t ysize = gridInqYsize(gridID);

  cdi_check_gridsize_int_limit("IEG", xsize * ysize);

  if (gridtype == GRID_GENERIC)
    {
      if ((ysize == 32 || ysize == 48 || ysize == 64 || ysize == 96 || ysize == 160) && (xsize == 2 * ysize || xsize == 1))
        {
          gridtype = GRID_GAUSSIAN;
          gridChangeType(gridID, gridtype);
        }
      else if ((xsize == 1 && ysize == 1) || (xsize == 0 && ysize == 0))
        {
          gridtype = GRID_LONLAT;
          gridChangeType(gridID, gridtype);
        }
      else if (gridInqXvals(gridID, NULL) && gridInqYvals(gridID, NULL))
        {
          gridtype = GRID_LONLAT;
          gridChangeType(gridID, gridtype);
        }
    }
  else if (gridtype == GRID_CURVILINEAR)
    {
      gridtype = GRID_LONLAT;
    }

  bool lrotated = (projtype == CDI_PROJ_RLL);

  if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || projtype == CDI_PROJ_RLL)
    {
      double xfirst = 0.0, xlast = 0.0, xinc = 0.0;
      double yfirst = 0.0, ylast = 0.0, yinc = 0.0;

      if (xsize == 0)
        xsize = 1;
      else
        {
          xfirst = gridInqXval(gridID, 0);
          xlast = gridInqXval(gridID, xsize - 1);
          xinc = gridInqXinc(gridID);
        }

      if (ysize == 0)
        ysize = 1;
      else
        {
          yfirst = gridInqYval(gridID, 0);
          ylast = gridInqYval(gridID, ysize - 1);
          yinc = gridInqYinc(gridID);
        }

      if (gridtype == GRID_GAUSSIAN)
        IEG_G_GridType(gdb) = 4;
      else if (lrotated)
        IEG_G_GridType(gdb) = 10;
      else
        IEG_G_GridType(gdb) = 0;

      double resfac = calc_resfac(xfirst, xlast, xinc, yfirst, ylast, yinc);
      int iresfac = (int) resfac;
      if (iresfac == 1000) iresfac = 0;

      IEG_G_ResFac(gdb) = iresfac;

      IEG_G_NumLon(gdb) = (int) xsize;
      IEG_G_NumLat(gdb) = (int) ysize;
      IEG_G_FirstLat(gdb) = (int) lround(yfirst * resfac);
      IEG_G_LastLat(gdb) = (int) lround(ylast * resfac);
      IEG_G_FirstLon(gdb) = (int) lround(xfirst * resfac);
      IEG_G_LastLon(gdb) = (int) lround(xlast * resfac);
      IEG_G_LonIncr(gdb) = (int) lround(xinc * resfac);
      if (fabs(xinc * resfac - IEG_G_LonIncr(gdb)) > FLT_EPSILON) IEG_G_LonIncr(gdb) = 0;

      if (gridtype == GRID_GAUSSIAN)
        IEG_G_LatIncr(gdb) = (int) ysize / 2;
      else
        {
          IEG_G_LatIncr(gdb) = (int) lround(yinc * resfac);
          if (fabs(yinc * resfac - IEG_G_LatIncr(gdb)) > FLT_EPSILON) IEG_G_LatIncr(gdb) = 0;

          if (IEG_G_LatIncr(gdb) < 0) IEG_G_LatIncr(gdb) = -IEG_G_LatIncr(gdb);
        }

      if (IEG_G_NumLon(gdb) > 1 && IEG_G_NumLat(gdb) == 1)
        if (IEG_G_LonIncr(gdb) != 0 && IEG_G_LatIncr(gdb) == 0) IEG_G_LatIncr(gdb) = IEG_G_LonIncr(gdb);

      if (IEG_G_NumLon(gdb) == 1 && IEG_G_NumLat(gdb) > 1)
        if (IEG_G_LonIncr(gdb) == 0 && IEG_G_LatIncr(gdb) != 0) IEG_G_LonIncr(gdb) = IEG_G_LatIncr(gdb);

      IEG_G_ResFlag(gdb) = (IEG_G_LatIncr(gdb) == 0 || IEG_G_LonIncr(gdb) == 0) ? 0 : 128;

      if (lrotated)
        {
          double xpole = 0.0, ypole = 0.0, angle = 0.0;
          gridInqParamRLL(gridID, &xpole, &ypole, &angle);

          IEG_G_LatSP(gdb) = -(int) lround(ypole * resfac);
          IEG_G_LonSP(gdb) = (int) lround((xpole + 180) * resfac);
          IEG_G_Size(gdb) = 42;
        }
      else
        {
          IEG_G_Size(gdb) = 32;
        }
    }
  else
    {
      Error("Unsupported grid type: %s", gridNamePtr(gridtype));
    }

  IEG_G_ScanFlag(gdb) = 64;
}

static void
pdbDefLevel(int *pdb, int leveltype, int level1, int level2)
{
  IEG_P_LevelType(pdb) = leveltype;
  IEG_P_Level1(pdb) = level1;
  IEG_P_Level2(pdb) = level2;
}

static void
iegDefLevel(int *pdb, int *gdb, double *vct, int zaxisID, int levelID)
{
  double level;

  int leveltype = zaxisInqType(zaxisID);
  if (leveltype == ZAXIS_GENERIC)
    {
      Warning("Changed zaxis type from %s to %s", zaxisNamePtr(leveltype), zaxisNamePtr(ZAXIS_PRESSURE));
      leveltype = ZAXIS_PRESSURE;
      zaxisChangeType(zaxisID, leveltype);
      cdiDefKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, "Pa");
    }

  // IEG_G_NumVCP(gdb) = 0;

  switch (leveltype)
    {
    case ZAXIS_SURFACE:
      {
        pdbDefLevel(pdb, IEG_LTYPE_SURFACE, 0, (int) (zaxisInqLevel(zaxisID, levelID)));
        break;
      }
    case ZAXIS_HYBRID:
      {
        if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
          pdbDefLevel(pdb, IEG_LTYPE_HYBRID_LAYER, (int) (zaxisInqLbound(zaxisID, levelID)),
                      (int) (zaxisInqUbound(zaxisID, levelID)));
        else
          pdbDefLevel(pdb, IEG_LTYPE_HYBRID, 0, (int) (zaxisInqLevel(zaxisID, levelID)));

        int vctsize = zaxisInqVctSize(zaxisID);
        if (vctsize > 100)
          {
            static bool vct_warning = true;
            // IEG_G_NumVCP(gdb) = 0;
            if (vct_warning)
              {
                Warning("VCT size of %d is too large (maximum is 100). Set to 0!", vctsize);
                vct_warning = false;
              }
          }
        else
          {
            IEG_G_Size(gdb) += (vctsize * 4);
            memcpy(vct, zaxisInqVctPtr(zaxisID), (size_t) vctsize / 2 * sizeof(double));
            memcpy(vct + 50, zaxisInqVctPtr(zaxisID) + vctsize / 2, (size_t) vctsize / 2 * sizeof(double));
          }
        break;
      }
    case ZAXIS_PRESSURE:
      {
        level = zaxisInqLevel(zaxisID, levelID);
        if (level < 0) Warning("pressure level of %f Pa is below 0.", level);

        char units[CDI_MAX_NAME];
        int length = CDI_MAX_NAME;
        cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_UNITS, units, &length);
        if (memcmp(units, "hPa", 3) == 0 || memcmp(units, "mb", 2) == 0) level = level * 100;

        double dum;
        int ilevel = (int) level;
        if (level < 32768 && (level < 100 || modf(level / 100, &dum) > 0))
          pdbDefLevel(pdb, IEG_LTYPE_99, 0, ilevel);
        else
          pdbDefLevel(pdb, IEG_LTYPE_ISOBARIC, 0, ilevel / 100);

        break;
      }
    case ZAXIS_HEIGHT:
      {
        level = zaxisInqLevel(zaxisID, levelID);
        pdbDefLevel(pdb, IEG_LTYPE_HEIGHT, 0, (int) level);
        break;
      }
    case ZAXIS_ALTITUDE:
      {
        level = zaxisInqLevel(zaxisID, levelID);
        pdbDefLevel(pdb, IEG_LTYPE_ALTITUDE, 0, (int) level);
        break;
      }
    case ZAXIS_DEPTH_BELOW_LAND:
      {
        if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
          pdbDefLevel(pdb, IEG_LTYPE_LANDDEPTH_LAYER, (int) (zaxisInqLbound(zaxisID, levelID)),
                      (int) (zaxisInqUbound(zaxisID, levelID)));
        else
          pdbDefLevel(pdb, IEG_LTYPE_LANDDEPTH, 0, (int) (zaxisInqLevel(zaxisID, levelID)));

        break;
      }
    case ZAXIS_DEPTH_BELOW_SEA:
      {
        level = zaxisInqLevel(zaxisID, levelID);
        pdbDefLevel(pdb, IEG_LTYPE_SEADEPTH, 0, (int) level);
        break;
      }
    case ZAXIS_ISENTROPIC:
      {
        level = zaxisInqLevel(zaxisID, levelID);
        pdbDefLevel(pdb, 113, 0, (int) level);
        break;
      }
    default:
      {
        Error("Unsupported zaxis type: %s", zaxisNamePtr(leveltype));
        break;
      }
    }
}

void
iegCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  streamFCopyRecord(streamptr2, streamptr1, "IEG");
}

void
iegDefRecord(stream_t *streamptr)
{
  Record *record = streamptr->record;

  int vlistID = streamptr->vlistID;
  int byteorder = streamptr->byteorder;

  int varID = record->varID;
  int levelID = record->levelID;
  int tsID = streamptr->curTsID;
  CdiDateTime vDateTime = streamptr->tsteps[tsID].taxis.vDateTime;

  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  iegrec_t *iegp = (iegrec_t *) record->objectp;
  iegInitMem(iegp);
  for (int i = 0; i < 37; i++) iegp->ipdb[i] = -1;

  iegp->byteswap = getByteswap(byteorder);

  int param = vlistInqVarParam(vlistID, varID);
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  IEG_P_Parameter(iegp->ipdb) = pnum;
  if (pdis == 255) IEG_P_CodeTable(iegp->ipdb) = pcat;

  iegDefTime(iegp->ipdb, vDateTime, vlistInqTaxis(vlistID));
  iegDefGrid(iegp->igdb, gridID);
  iegDefLevel(iegp->ipdb, iegp->igdb, iegp->vct, zaxisID, levelID);

  iegp->dprec = iegDefDatatype(record->prec);
}

static double
calc_refvalSP(size_t size, const float *data)
{
  float refval = data[0];
  for (size_t i = 1; i < size; ++i)
    if (data[i] < refval) refval = data[i];
  return refval;
}

static double
calc_refvalDP(size_t size, const double *data)
{
  double refval = data[0];
  for (size_t i = 1; i < size; ++i)
    if (data[i] < refval) refval = data[i];
  return refval;
}

static void
ieg_write_recordSP(stream_t *streamptr, const float *data)
{
  Record *record = streamptr->record;
  iegrec_t *iegp = (iegrec_t *) record->objectp;

  int fileID = streamptr->fileID;
  size_t gridsize = gridInqSize(record->gridID);

  iegp->refval = calc_refvalSP(gridsize, data);

  iegDefDataSP(iegp, data);
  iegWrite(fileID, iegp);
}

static void
ieg_write_recordDP(stream_t *streamptr, const double *data)
{
  Record *record = streamptr->record;
  iegrec_t *iegp = (iegrec_t *) record->objectp;

  int fileID = streamptr->fileID;
  size_t gridsize = gridInqSize(record->gridID);

  iegp->refval = calc_refvalDP(gridsize, data);

  iegDefDataDP(iegp, data);
  iegWrite(fileID, iegp);
}

void
ieg_write_record(stream_t *streamptr, int memtype, const void *data)
{
  if (memtype == MEMTYPE_DOUBLE)
    ieg_write_recordDP(streamptr, (const double *) data);
  else
    ieg_write_recordSP(streamptr, (const float *) data);
}

static void
iegAddRecord(stream_t *streamptr, int param, int *pdb, int *gdb, double *vct, size_t recsize, off_t position, int prec)
{
  int vlistID = streamptr->vlistID;
  int tsID = streamptr->curTsID;
  int recID = recordNewEntry(streamptr, tsID);
  record_t *record = &streamptr->tsteps[tsID].records[recID];

  int level1, level2;
  if (IEG_P_LevelType(pdb) == IEG_LTYPE_HYBRID_LAYER)
    {
      level1 = IEG_P_Level1(pdb);
      level2 = IEG_P_Level2(pdb);
    }
  else
    {
      level1 = IEG_P_Level2(pdb);
      level2 = 0;
      if (IEG_P_LevelType(pdb) == 100) level1 *= 100;
    }

  record->size = recsize;
  record->position = position;
  record->param = param;
  record->ilevel = level1;
  record->ilevel2 = level2;
  record->ltype = IEG_P_LevelType(pdb);

  int gridtype = (IEG_G_GridType(gdb) == 0)    ? GRID_LONLAT
                 : (IEG_G_GridType(gdb) == 10) ? GRID_PROJECTION
                 : (IEG_G_GridType(gdb) == 4)  ? GRID_GAUSSIAN
                                               : GRID_GENERIC;

  grid_t *grid = (grid_t *) Malloc(sizeof(*grid));
  grid_init(grid);
  cdiGridTypeInit(grid, gridtype, IEG_G_NumLon(gdb) * IEG_G_NumLat(gdb));
  size_t xsize = (size_t) IEG_G_NumLon(gdb);
  size_t ysize = (size_t) IEG_G_NumLat(gdb);
  grid->x.size = xsize;
  grid->y.size = ysize;
  grid->x.inc = 0;
  grid->y.inc = 0;
  grid->x.flag = 0;

  int iresfac = IEG_G_ResFac(gdb);
  if (iresfac == 0) iresfac = 1000;
  double resfac = 1. / (double) iresfac;

  // if ( IEG_G_FirstLon != 0 || IEG_G_LastLon != 0 )
  {
    if (xsize > 1)
      {
        if (IEG_G_ResFlag(gdb) && IEG_G_LonIncr(gdb) > 0)
          grid->x.inc = IEG_G_LonIncr(gdb) * resfac;
        else
          grid->x.inc = (IEG_G_LastLon(gdb) - IEG_G_FirstLon(gdb)) * resfac / (xsize - 1);

        // correct xinc if necessary
        if (IEG_G_FirstLon(gdb) == 0 && IEG_G_LastLon(gdb) > 354000)
          {
            double xinc = 360. / xsize;
            // FIXME: why not use grid->x.inc != xinc as condition?
            if (fabs(grid->x.inc - xinc) > 0.0)
              {
                grid->x.inc = xinc;
                if (CDI_Debug) Message("set xinc to %g", grid->x.inc);
              }
          }
      }
    grid->x.first = IEG_G_FirstLon(gdb) * resfac;
    grid->x.last = IEG_G_LastLon(gdb) * resfac;
    grid->x.flag = 2;
  }
  grid->y.flag = 0;
  // if ( IEG_G_FirstLat != 0 || IEG_G_LastLat != 0 )
  {
    if (ysize > 1)
      {
        if (IEG_G_ResFlag(gdb) && IEG_G_LatIncr(gdb) > 0)
          grid->y.inc = IEG_G_LatIncr(gdb) * resfac;
        else
          grid->y.inc = (IEG_G_LastLat(gdb) - IEG_G_FirstLat(gdb)) * resfac / (ysize - 1);
      }
    grid->y.first = IEG_G_FirstLat(gdb) * resfac;
    grid->y.last = IEG_G_LastLat(gdb) * resfac;
    grid->y.flag = 2;
  }

  double xpole = 0, ypole = 0;
  if (IEG_G_GridType(gdb) == 10)
    {
      xpole = IEG_G_LonSP(gdb) * resfac - 180;
      ypole = -IEG_G_LatSP(gdb) * resfac;
      grid->projtype = CDI_PROJ_RLL;
    }

  struct addIfNewRes gridAdded = cdiVlistAddGridIfNew(vlistID, grid, 0);
  int gridID = gridAdded.Id;
  if (!gridAdded.isNew)
    {
      grid_free(grid);
      Free(grid);
    }
  else if (gridtype == GRID_PROJECTION)
    gridDefParamRLL(gridID, xpole, ypole, 0);

  int leveltype = iegGetZaxisType(IEG_P_LevelType(pdb));
  if (leveltype == ZAXIS_HYBRID)
    {
      double tmpvct[100];
      size_t vctsize = (size_t) IEG_G_NumVCP(gdb);

      for (size_t i = 0; i < vctsize / 2; i++) tmpvct[i] = vct[i];
      for (size_t i = 0; i < vctsize / 2; i++) tmpvct[i + vctsize / 2] = vct[i + 50];

      varDefVCT(vctsize, tmpvct);
    }

  int lbounds = IEG_P_LevelType(pdb) == IEG_LTYPE_HYBRID_LAYER ? 1 : 0;

  int datatype = iegInqDatatype(prec);

  int varID, levelID = 0;
  varAddRecord(recID, param, gridID, leveltype, lbounds, level1, level2, 0, 0, datatype, &varID, &levelID, TSTEP_INSTANT, 0, -1,
               NULL, NULL, NULL, NULL);

  record->varID = (short) varID;
  record->levelID = levelID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;

  if (CDI_Debug) Message("varID = %d gridID = %d levelID = %d", varID, gridID, levelID);
}

#if 0
static
void iegCmpRecord(stream_t *streamptr, int tsID, int recID, off_t position, int param,
		  int level, size_t xsize, size_t ysize)
{
  int varID = 0;
  int levelID = 0;

  record_t *record  = &streamptr->tsteps[tsID].records[recID];

  if (param != (*record).param || level != (*record).ilevel) Error("inconsistent timestep");

  (*record).position = position;
  /*
  varID   = (*record).varID;
  levelID = (*record).levelID;

  streamptr->vars[varID].level[levelID] = recID;

  streamptr->tsteps[tsID].nallrecs++;
  streamptr->nrecs++;
  */
  if (CDI_Debug) Message("varID = %d levelID = %d", varID, levelID);
}
#endif

static CdiDateTime
iegDateTime(const int *pdb)
{
  int ryear = IEG_P_Year(pdb);
  int rmonth = IEG_P_Month(pdb);
  int rday = IEG_P_Day(pdb);

  int rhour = IEG_P_Hour(pdb);
  int rminute = IEG_P_Minute(pdb);

  if (rminute == -1) rminute = 0;

  CdiDateTime cdiDateTime;
  cdiDateTime.date = cdiDate_encode(ryear, rmonth, rday);
  cdiDateTime.time = cdiTime_encode(rhour, rminute, 0, 0);
  return cdiDateTime;
}

static void
iegScanTimestep1(stream_t *streamptr)
{
  CdiDateTime datetime0;
  cdiDateTime_init(&datetime0);
  off_t recpos;
  iegrec_t *iegp = (iegrec_t *) streamptr->record->objectp;

  streamptr->curTsID = 0;

  int tsID = tstepsNewEntry(streamptr);
  if (tsID != 0) Error("Internal problem! tstepsNewEntry returns %d", tsID);
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  int fileID = streamptr->fileID;

  int nrecs = 0;
  while (true)
    {
      recpos = fileGetPos(fileID);
      if (iegRead(fileID, iegp) != 0)
        {
          streamptr->ntsteps = 1;
          break;
        }

      size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

      int prec = iegp->dprec;
      int rcode = IEG_P_Parameter(iegp->ipdb);
      int tabnum = IEG_P_CodeTable(iegp->ipdb);
      int param = cdiEncodeParam(rcode, tabnum, 255);
      int rlevel = (IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER) ? IEG_P_Level1(iegp->ipdb) : IEG_P_Level2(iegp->ipdb);

      if (IEG_P_LevelType(iegp->ipdb) == 100) rlevel *= 100;

      CdiDateTime datetime = iegDateTime(iegp->ipdb);

      if (nrecs == 0)
        {
          datetime0 = datetime;
          taxis->vDateTime = datetime;
        }
      else
        {
          record_t *records = streamptr->tsteps[tsID].records;
          for (int recID = 0; recID < nrecs; recID++)
            if (param == records[recID].param && rlevel == records[recID].ilevel) goto tstepScanLoopFinished;

          if (cdiDateTime_isNE(datetime, datetime0)) Warning("Inconsistent verification time for param %d level %d", param, rlevel);
        }

      nrecs++;

      if (CDI_Debug) Message("%4d%8d%4d%8d %s", nrecs, (int) recpos, param, rlevel, CdiDateTime_string(datetime));

      iegAddRecord(streamptr, param, iegp->ipdb, iegp->igdb, iegp->vct, recsize, recpos, prec);
    }

tstepScanLoopFinished:
  streamptr->rtsteps = 1;

  cdi_generate_vars(streamptr);

  int taxisID = taxisCreate(TAXIS_ABSOLUTE);
  taxis->type = TAXIS_ABSOLUTE;
  taxis->rDateTime = taxis->vDateTime;

  int vlistID = streamptr->vlistID;
  vlistDefTaxis(vlistID, taxisID);

  vlist_check_contents(vlistID);

  streamScanResizeRecords1(streamptr);

  streamScanTsFixNtsteps(streamptr, recpos);
  streamScanTimeConstAdjust(streamptr, taxis);
}

static int
iegScanTimestep2(stream_t *streamptr)
{
  off_t recpos = 0;
  iegrec_t *iegp = (iegrec_t *) streamptr->record->objectp;

  streamptr->curTsID = 1;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  int tsID = streamptr->rtsteps;
  if (tsID != 1) Error("Internal problem! unexpected timestep %d", tsID + 1);

  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

  cdi_create_records(streamptr, tsID);
  record_t *records = streamptr->tsteps[tsID].records;

  int nrecords = streamScanInitRecords2(streamptr);

  for (int rindex = 0; rindex <= nrecords; rindex++)
    {
      recpos = fileGetPos(fileID);
      if (iegRead(fileID, iegp) != 0)
        {
          streamptr->ntsteps = 2;
          break;
        }

      size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

      int rcode = IEG_P_Parameter(iegp->ipdb);
      int tabnum = IEG_P_CodeTable(iegp->ipdb);
      int param = cdiEncodeParam(rcode, tabnum, 255);
      int rlevel = (IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER) ? IEG_P_Level1(iegp->ipdb) : IEG_P_Level2(iegp->ipdb);

      if (IEG_P_LevelType(iegp->ipdb) == 100) rlevel *= 100;

      CdiDateTime datetime = iegDateTime(iegp->ipdb);

      if (rindex == 0)
        {
          taxis->type = TAXIS_ABSOLUTE;
          taxis->vDateTime = datetime;
        }

      bool nextstep = false;
      int recID = 0;
      for (recID = 0; recID < nrecords; recID++)
        {
          if (param == records[recID].param && rlevel == records[recID].ilevel)
            {
              if (records[recID].used)
                {
                  nextstep = true;
                }
              else
                {
                  records[recID].used = true;
                  streamptr->tsteps[tsID].recIDs[rindex] = recID;
                }
              break;
            }
        }
      if (recID == nrecords)
        {
          char paramstr[32];
          cdiParamToString(param, paramstr, sizeof(paramstr));
          Warning("param %s level %d not defined at timestep 1", paramstr, rlevel);
          return CDI_EUFSTRUCT;
        }

      if (nextstep) break;

      if (CDI_Debug) Message("%4d%8d%4d%8d %s", rindex + 1, (int) recpos, param, rlevel, CdiDateTime_string(datetime));

      if (param != records[recID].param || rlevel != records[recID].ilevel)
        {
          Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                  records[recID].ilevel, rlevel);
          return CDI_EUFSTRUCT;
        }

      records[recID].position = recpos;
      records[recID].size = recsize;
    }

  int nrecs = 0;
  for (int recID = 0; recID < nrecords; recID++)
    {
      if (records[recID].used)
        nrecs++;
      else
        vlistDefVarTimetype(vlistID, records[recID].varID, TIME_CONSTANT);
    }
  streamptr->tsteps[tsID].nrecs = nrecs;

  streamptr->rtsteps = 2;

  streamScanTsFixNtsteps(streamptr, recpos);

  return 0;
}

int
iegInqContents(stream_t *streamptr)
{
  streamptr->curTsID = 0;

  iegScanTimestep1(streamptr);

  int status = (streamptr->ntsteps == -1) ? iegScanTimestep2(streamptr) : 0;

  fileSetPos(streamptr->fileID, 0, SEEK_SET);

  return status;
}

static long
iegScanTimestep(stream_t *streamptr)
{
  off_t recpos = 0;
  iegrec_t *iegp = (iegrec_t *) streamptr->record->objectp;

  if (streamptr->rtsteps == 0) Error("Internal problem! Missing contents.");

  int tsID = streamptr->rtsteps;
  taxis_t *taxis = &streamptr->tsteps[tsID].taxis;

  int nrecs = 0;
  if (streamptr->tsteps[tsID].recordSize == 0)
    {
      cdi_create_records(streamptr, tsID);
      record_t *records = streamptr->tsteps[tsID].records;

      nrecs = streamScanInitRecords(streamptr, tsID);

      int fileID = streamptr->fileID;

      fileSetPos(fileID, streamptr->tsteps[tsID].position, SEEK_SET);

      for (int rindex = 0; rindex <= nrecs; rindex++)
        {
          recpos = fileGetPos(fileID);
          if (iegRead(fileID, iegp) != 0)
            {
              streamptr->ntsteps = streamptr->rtsteps + 1;
              break;
            }

          size_t recsize = (size_t) (fileGetPos(fileID) - recpos);

          int rcode = IEG_P_Parameter(iegp->ipdb);
          int tabnum = IEG_P_CodeTable(iegp->ipdb);
          int param = cdiEncodeParam(rcode, tabnum, 255);
          int rlevel
              = (IEG_P_LevelType(iegp->ipdb) == IEG_LTYPE_HYBRID_LAYER) ? IEG_P_Level1(iegp->ipdb) : IEG_P_Level2(iegp->ipdb);

          if (IEG_P_LevelType(iegp->ipdb) == 100) rlevel *= 100;

          CdiDateTime datetime = iegDateTime(iegp->ipdb);

          // if ( rindex == nrecs ) break; gcc-4.5 internal compiler error
          if (rindex == nrecs) continue;
          int recID = streamptr->tsteps[tsID].recIDs[rindex];

          if (rindex == 0)
            {
              taxis->type = TAXIS_ABSOLUTE;
              taxis->vDateTime = datetime;
            }

          if (param != records[recID].param || rlevel != records[recID].ilevel)
            {
              Message("tsID = %d recID = %d param = %3d new %3d  level = %3d new %3d", tsID, recID, records[recID].param, param,
                      records[recID].ilevel, rlevel);
              Error("Invalid, unsupported or inconsistent record structure");
            }

          records[recID].position = recpos;
          records[recID].size = recsize;

          if (CDI_Debug) Message("%4d%8d%4d%8d %s", rindex, (int) recpos, param, rlevel, CdiDateTime_string(datetime));
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
    }

  if (nrecs > 0 && nrecs < streamptr->tsteps[tsID].nrecs)
    {
      Warning("Incomplete timestep. Stop scanning at timestep %d.", tsID);
      streamptr->ntsteps = tsID;
    }

  return streamptr->ntsteps;
}

int
iegInqTimestep(stream_t *streamptr, int tsID)
{
  if (tsID == 0 && streamptr->rtsteps == 0) Error("Call to cdiInqContents missing!");

  if (CDI_Debug) Message("tsID = %d rtsteps = %d", tsID, streamptr->rtsteps);

  long ntsteps = CDI_UNDEFID;
  while ((tsID + 1) > streamptr->rtsteps && ntsteps == CDI_UNDEFID) ntsteps = iegScanTimestep(streamptr);

  int nrecs = 0;
  if (!(tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID))
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
    }

  return nrecs;
}

void
iegReadVarSliceDP(stream_t *streamptr, int varID, int levID, double *data, size_t *numMissVals)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  void *iegp = streamptr->record->objectp;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  double missval = vlistInqVarMissval(vlistID, varID);
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  int tsid = streamptr->curTsID;

  off_t currentfilepos = fileGetPos(fileID);

  int recID = streamptr->vars[varID].recordTable[0].recordID[levID];
  off_t recpos = streamptr->tsteps[tsid].records[recID].position;
  fileSetPos(fileID, recpos, SEEK_SET);
  iegRead(fileID, iegp);
  iegInqDataDP(iegp, data);

  fileSetPos(fileID, currentfilepos, SEEK_SET);

  *numMissVals = get_num_missvalsDP(gridsize, data, missval);
}

void
iegReadVarDP(stream_t *streamptr, int varID, double *data, size_t *numMissVals)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevs = (size_t) streamptr->vars[varID].recordTable[0].nlevs;

  for (size_t levID = 0; levID < nlevs; levID++)
    iegReadVarSliceDP(streamptr, varID, (int) levID, &data[levID * gridsize], numMissVals);
}

void
iegWriteVarSliceDP(stream_t *streamptr, int varID, int levID, const double *data)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d  levID = %d", streamptr->self, varID, levID);

  iegrec_t *iegp = (iegrec_t *) streamptr->record->objectp;
  iegInitMem(iegp);
  for (int i = 0; i < 37; i++) iegp->ipdb[i] = -1;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int tsID = streamptr->curTsID;
  CdiDateTime vDateTime = streamptr->tsteps[tsID].taxis.vDateTime;
  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  int param = vlistInqVarParam(vlistID, varID);
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  IEG_P_Parameter(iegp->ipdb) = pnum;
  if (pdis == 255) IEG_P_CodeTable(iegp->ipdb) = pcat;

  iegDefTime(iegp->ipdb, vDateTime, vlistInqTaxis(vlistID));
  iegDefGrid(iegp->igdb, gridID);
  iegDefLevel(iegp->ipdb, iegp->igdb, iegp->vct, zaxisID, levID);

  iegp->dprec = iegDefDatatype(vlistInqVarDatatype(vlistID, varID));

  size_t gridsize = gridInqSize(gridID);
  iegp->refval = calc_refvalDP(gridsize, data);

  iegDefDataDP(iegp, data);
  iegWrite(fileID, iegp);
}

void
iegWriteVarDP(stream_t *streamptr, int varID, const double *data)
{
  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
  size_t nlevels = (size_t) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));

  for (size_t levID = 0; levID < nlevels; levID++) iegWriteVarSliceDP(streamptr, varID, (int) levID, &data[levID * gridsize]);
}

#endif /* HAVE_LIBIEG */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
