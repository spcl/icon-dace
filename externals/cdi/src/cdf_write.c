#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBNETCDF

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "stream_cdf.h"
#include "cdf.h"
#include "cdf_int.h"
#include "vlist.h"
#include "vlist_var.h"

// #include <netcdf_filter.h>

void
cdf_def_var_filter(int ncid, int ncvarID, unsigned int id, size_t nparams, const unsigned int *params)
{
#ifdef HAVE_NETCDF4
  int status;
  if ((status = nc_def_var_filter(ncid, ncvarID, id, nparams, params)))
    {
      Message("filterId=%u  numParams=%zu", id, nparams);
      Error("nc_def_var_filter failed; %s", nc_strerror(status));
    }
#else
  static bool lwarn = true;
  if (lwarn)
    {
      lwarn = false;
      Warning("filter failed, NetCDF4 not available!");
    }
#endif
}

void
cdfDefVarDeflate(int ncid, int ncvarID, int compLevel)
{
#ifdef HAVE_NETCDF4
  // Set chunking, shuffle, and deflate.
  int shuffle = (CDI_Shuffle > 0), deflate = 1;

  if (compLevel < 1 || compLevel > 9) compLevel = 1;

  int status;
  if ((status = nc_def_var_deflate(ncid, ncvarID, shuffle, deflate, compLevel)))
    {
      Error("nc_def_var_deflate failed; %s", nc_strerror(status));
    }
#else
  static bool lwarn = true;
  if (lwarn)
    {
      lwarn = false;
      Warning("Deflate compression failed, NetCDF4 not available!");
    }
#endif
}

void
cdfDefVarSzip(int ncid, int ncvarID, int pixels_per_block)
{
#ifdef HAVE_NC_DEF_VAR_SZIP
  // Set options_mask.
  /*
    H5_SZIP_ALLOW_K13_OPTION_MASK   1
    H5_SZIP_CHIP_OPTION_MASK        2
    H5_SZIP_EC_OPTION_MASK          4
    H5_SZIP_NN_OPTION_MASK          32
    H5_SZIP_ALL_MASKS (H5_SZIP_CHIP_OPTION_MASK|H5_SZIP_EC_OPTION_MASK|H5_SZIP_NN_OPTION_MASK)
  */
  int options_mask = 38;
  int status;
  if ((status = nc_def_var_szip(ncid, ncvarID, options_mask, pixels_per_block)))
    Error("nc_def_var_szip failed; %s", nc_strerror(status));
#else
  (void) ncid;
  (void) ncvarID;
  (void) pixels_per_block;
  static bool lwarn = true;
  if (lwarn)
    {
      lwarn = false;
      Warning("Szip compression failed, NetCDF4/szlib not available!");
    }
#endif
}

#ifdef HAVE_NETCDF4
static nc_type
cdfTypeComplexFloat(stream_t *streamptr)
{
  if (streamptr->nc_complex_float_id == CDI_UNDEFID)
    {
      typedef struct complex_float
      {
        float r, i;
      } complex_float;
      int fileID = streamptr->fileID;
      int nc_complex_id;
      int status;
      status = nc_def_compound(fileID, sizeof(complex_float), "complex_float", &nc_complex_id);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      status = nc_insert_compound(fileID, nc_complex_id, "r", NC_COMPOUND_OFFSET(complex_float, r), NC_FLOAT);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      status = nc_insert_compound(fileID, nc_complex_id, "i", NC_COMPOUND_OFFSET(complex_float, i), NC_FLOAT);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      streamptr->nc_complex_float_id = nc_complex_id;
    }

  return (nc_type) streamptr->nc_complex_float_id;
}

static nc_type
cdfTypeComplexDouble(stream_t *streamptr)
{
  if (streamptr->nc_complex_double_id == CDI_UNDEFID)
    {
      typedef struct complex_double
      {
        double r, i;
      } complex_double;
      int fileID = streamptr->fileID;
      int nc_complex_id;
      int status;
      status = nc_def_compound(fileID, sizeof(complex_double), "complex_double", &nc_complex_id);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      status = nc_insert_compound(fileID, nc_complex_id, "r", NC_COMPOUND_OFFSET(complex_double, r), NC_DOUBLE);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      status = nc_insert_compound(fileID, nc_complex_id, "i", NC_COMPOUND_OFFSET(complex_double, i), NC_DOUBLE);
      if (status != NC_NOERR) Error("%s", nc_strerror(status));
      streamptr->nc_complex_double_id = nc_complex_id;
    }

  return (nc_type) streamptr->nc_complex_double_id;
}
#endif

nc_type
cdfDefDatatype(int datatype, stream_t *streamptr)
{
  nc_type xtype = NC_FLOAT;

  // clang-format off
  if (streamptr->filetype == CDI_FILETYPE_NC4)
    {
      if      (datatype == CDI_DATATYPE_INT8  ) xtype = NC_BYTE;
      else if (datatype == CDI_DATATYPE_INT16 ) xtype = NC_SHORT;
      else if (datatype == CDI_DATATYPE_INT32 ) xtype = NC_INT;
#ifdef  HAVE_NETCDF4
      else if (datatype == CDI_DATATYPE_UINT8 ) xtype = NC_UBYTE;
      else if (datatype == CDI_DATATYPE_UINT16) xtype = NC_USHORT;
      else if (datatype == CDI_DATATYPE_UINT32) xtype = NC_UINT;
      else if (datatype == CDI_DATATYPE_CPX32 ) xtype = cdfTypeComplexFloat(streamptr);
      else if (datatype == CDI_DATATYPE_CPX64 ) xtype = cdfTypeComplexDouble(streamptr);
#else
      else if (datatype == CDI_DATATYPE_UINT8 ) xtype = NC_SHORT;
      else if (datatype == CDI_DATATYPE_UINT16) xtype = NC_INT;
      else if (datatype == CDI_DATATYPE_UINT32) xtype = NC_INT;
      else if (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64)
        Error("CDI library does not support complex numbers with NetCDF4 classic!");
#endif
      else if (datatype == CDI_DATATYPE_FLT64 ) xtype = NC_DOUBLE;
      else if (datatype == CDI_DATATYPE_FLT32 ) xtype = NC_FLOAT;
    }
  else
    {
      if      (datatype == CDI_DATATYPE_INT8  ) xtype = NC_BYTE;
      else if (datatype == CDI_DATATYPE_INT16 ) xtype = NC_SHORT;
      else if (datatype == CDI_DATATYPE_INT32 ) xtype = NC_INT;
      else if (datatype == CDI_DATATYPE_UINT8 ) xtype = NC_SHORT;
      else if (datatype == CDI_DATATYPE_UINT16) xtype = NC_INT;
      else if (datatype == CDI_DATATYPE_UINT32) xtype = NC_INT;
      else if (datatype == CDI_DATATYPE_FLT64 ) xtype = NC_DOUBLE;
      else if (datatype == CDI_DATATYPE_FLT32 ) xtype = NC_FLOAT;
      else if (datatype == CDI_DATATYPE_CPX32 || datatype == CDI_DATATYPE_CPX64)
        Error("CDI library does not support complex numbers with NetCDF classic!");
    }
  // clang-format on

  return xtype;
}

static void
cdfDefVarMissval(stream_t *streamptr, int varID, int dtype, int lcheck)
{
  if (streamptr->vars[varID].defmiss == false)
    {
      int vlistID = streamptr->vlistID;
      int fileID = streamptr->fileID;
      int ncvarID = streamptr->vars[varID].ncvarid;
      double missval = vlistInqVarMissval(vlistID, varID);

      if (lcheck && streamptr->ncmode == 2) cdf_redef(fileID);

      nc_type xtype = cdfDefDatatype(dtype, streamptr);
      if (xtype == NC_BYTE && missval > 127 && missval < 256) xtype = NC_INT;

      if (lcheck == 0 || streamptr->ncmode != 2 || streamptr->filetype == CDI_FILETYPE_NC || streamptr->filetype == CDI_FILETYPE_NC2
          || streamptr->filetype == CDI_FILETYPE_NC5)
        cdf_put_att_double(fileID, ncvarID, "_FillValue", xtype, 1, &missval);

      cdf_put_att_double(fileID, ncvarID, "missing_value", xtype, 1, &missval);

      if (lcheck && streamptr->ncmode == 2) cdf_enddef(fileID, streamptr->self);

      streamptr->vars[varID].defmiss = true;
    }
}

static void
cdfDefInstituteGlobal(const stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int instID = vlistInqInstitut(vlistID);

  if (instID != CDI_UNDEFID)
    {
      const char *longname = institutInqLongnamePtr(instID);
      if (longname)
        {
          size_t len = strlen(longname);
          if (len > 0)
            {
              if (streamptr->ncmode == 2) cdf_redef(fileID);
              cdf_put_att_text(fileID, NC_GLOBAL, "institution", len, longname);
              if (streamptr->ncmode == 2) cdf_enddef(fileID, streamptr->self);
            }
        }
    }
}

static void
cdfDefSourceGlobal(const stream_t *streamptr)
{
  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;
  int modelID = vlistInqModel(vlistID);

  if (modelID != CDI_UNDEFID)
    {
      const char *longname = modelInqNamePtr(modelID);
      if (longname)
        {
          size_t len = strlen(longname);
          if (len > 0)
            {
              if (streamptr->ncmode == 2) cdf_redef(fileID);
              cdf_put_att_text(fileID, NC_GLOBAL, "source", len, longname);
              if (streamptr->ncmode == 2) cdf_enddef(fileID, streamptr->self);
            }
        }
    }
}

static inline void *
resizeBuf(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static void
cdfDefineCellMethods(stream_t *streamptr, int cdiID, int varID, int fileID, int ncvarID)
{
  taxis_t *taxis = &streamptr->tsteps[0].taxis;
  if (!taxis->hasBounds) return;

  int timeVarId = streamptr->basetime.ncvarid;
  char timeVarName[CDI_MAX_NAME];
  cdf_inq_varname(fileID, timeVarId, timeVarName);

  int stepType = vlistInqVarTsteptype(cdiID, varID);

  const char *cellMethod = NULL;
  // clang-format off
  if      (stepType == TSTEP_AVG)   cellMethod = "mean";
  else if (stepType == TSTEP_SUM)   cellMethod = "sum";
  else if (stepType == TSTEP_RANGE) cellMethod = "range";
  else if (stepType == TSTEP_MIN)   cellMethod = "minimum";
  else if (stepType == TSTEP_MAX)   cellMethod = "maximum";
  // clang-format on

  if (cellMethod)
    {
      const char *attname = "cell_methods";
      char atttxt[CDI_MAX_NAME + 10];
      snprintf(atttxt, sizeof(atttxt), "%s: %s", timeVarName, cellMethod);
      cdf_put_att_text(fileID, ncvarID, attname, strlen(atttxt), atttxt);
    }
}

static nc_type
int_datatype_to_xtype(int filetype, int datatype)
{
  // clang-format off
  if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C || filetype == CDI_FILETYPE_NCZARR)
    {
      return (datatype == CDI_DATATYPE_INT8)   ? NC_BYTE :
             (datatype == CDI_DATATYPE_INT16)  ? NC_SHORT :
#ifdef  HAVE_NETCDF4
             (datatype == CDI_DATATYPE_UINT8)  ? NC_UBYTE :
             (datatype == CDI_DATATYPE_UINT16) ? NC_USHORT :
             (datatype == CDI_DATATYPE_UINT32) ? NC_UINT :
#endif
                                                 NC_INT;
    }

  return (datatype == CDI_DATATYPE_INT8)   ? NC_BYTE :
         (datatype == CDI_DATATYPE_INT16)  ? NC_SHORT :
                                             NC_INT;
  // clang-format on
}

void
cdfDefineAttributes(int filetype, int cdiID, int varID, int fileID, int ncvarID)
{
  int atttype, attlen;
  char attname[CDI_MAX_NAME + 1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for (int iatt = 0; iatt < natts; ++iatt)
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      // if (attlen == 0) continue;

      if (atttype == CDI_DATATYPE_TXT)
        {
          size_t attSize = (size_t) attlen * sizeof(char);
          char *atttxt = (char *) resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          size_t len = (size_t) attlen;
          cdf_put_att_text(fileID, ncvarID, attname, len, atttxt);
        }
      else if (atttype == CDI_DATATYPE_INT8 || atttype == CDI_DATATYPE_UINT8 || atttype == CDI_DATATYPE_INT16
               || atttype == CDI_DATATYPE_UINT16 || atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32)
        {
          if (attlen == 0) continue;
          size_t attSize = (size_t) attlen * sizeof(int);
          int *attint = (int *) resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          size_t len = (size_t) attlen;
          cdf_put_att_int(fileID, ncvarID, attname, int_datatype_to_xtype(filetype, atttype), len, attint);
        }
      else if (atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64)
        {
          if (attlen == 0) continue;
          size_t attSize = (size_t) attlen * sizeof(double);
          double *attflt = (double *) resizeBuf(&attBuf, &attBufSize, attSize);
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          size_t len = (size_t) attlen;
          if (atttype == CDI_DATATYPE_FLT32)
            {
              float attflt_sp[8];
              float *pattflt_sp = (len > 8) ? (float *) malloc(len * sizeof(float)) : attflt_sp;
              for (size_t i = 0; i < len; ++i) pattflt_sp[i] = (float) attflt[i];
              cdf_put_att_float(fileID, ncvarID, attname, NC_FLOAT, len, pattflt_sp);
              if (len > 8) free(pattflt_sp);
            }
          else
            cdf_put_att_double(fileID, ncvarID, attname, NC_DOUBLE, len, attflt);
        }
    }

  if (attBuf) Free(attBuf);
}

static void
cdfDefineInstituteName(int vlistID, int varID, int fileID, int ncvarID)
{
  int instID = vlistInqVarInstitut(vlistID, varID);
  if (instID != CDI_UNDEFID)
    {
      const char *name = institutInqNamePtr(instID);
      if (name) cdf_put_att_text(fileID, ncvarID, "institution", strlen(name), name);
    }
}

static void
cdfDefGlobalAtts(stream_t *streamptr)
{
  if (streamptr->globalatts) return;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  cdfDefSourceGlobal(streamptr);
  cdfDefInstituteGlobal(streamptr);

  int natts;
  cdiInqNatts(vlistID, CDI_GLOBAL, &natts);

  if (natts > 0 && streamptr->ncmode == 2) cdf_redef(fileID);

  cdfDefineAttributes(streamptr->filetype, vlistID, CDI_GLOBAL, fileID, NC_GLOBAL);

  if (natts > 0 && streamptr->ncmode == 2) cdf_enddef(fileID, streamptr->self);

  streamptr->globalatts = 1;
}

static void
cdf_get_gmapvarname(int gridID, char *gmapvarname)
{
  int length = CDI_MAX_NAME;
  int pgridID = gridID;
  cdiInqKeyString(pgridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, gmapvarname, &length);

  if (!gmapvarname[0])
    {
      length = CDI_MAX_NAME;
      pgridID = gridInqProj(gridID);
      if (pgridID != CDI_UNDEFID) cdiInqKeyString(pgridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, gmapvarname, &length);
    }
}

static int
nc_grid_index(stream_t *streamptr, int gridID)
{
  int index = 0;
  int vlistID = streamptr->vlistID;
  int ngrids = vlistNgrids(vlistID);
  for (index = 0; index < ngrids; ++index)
    if (streamptr->ncgrid[index].gridID == gridID) break;

  assert(index < ngrids);

  return index;
}

// convert NetCDF xtype to pixels_per_block
static int
xtype2ppb(nc_type xtype)
{
  int ppb = 32;

  // clang-format off
  if      (xtype == NC_BYTE)   ppb = 8;
  else if (xtype == NC_SHORT)  ppb = 16;
#ifdef  HAVE_NETCDF4
  else if (xtype == NC_UBYTE)  ppb = 8;
  else if (xtype == NC_USHORT) ppb = 16;
#endif
  // clang-format on

  return ppb;
}

static void
cdfDefVarFilter(const stream_t *s, int ncvarID)
{
  if (s->filterId != 0)
    {
      if (s->filetype == CDI_FILETYPE_NC4 || s->filetype == CDI_FILETYPE_NC4C || s->filetype == CDI_FILETYPE_NCZARR)
        {
          cdf_def_var_filter(s->fileID, ncvarID, s->filterId, s->numParams, s->params);
        }
      else
        {
          static bool lwarn = true;
          if (lwarn)
            {
              lwarn = false;
              Warning("Filter is only available for NetCDF4!");
            }
        }
    }
}

static void
cdfDefVarCompression(const stream_t *streamptr, int ncvarID, nc_type xtype)
{
  if (streamptr->comptype == CDI_COMPRESS_ZIP)
    {
      if (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C
          || streamptr->filetype == CDI_FILETYPE_NCZARR)
        {
          cdfDefVarDeflate(streamptr->fileID, ncvarID, streamptr->complevel);
        }
      else
        {
          static bool lwarn = true;
          if (lwarn)
            {
              lwarn = false;
              Warning("Deflate compression is only available for NetCDF4!");
            }
        }
    }
  /*
  else if (streamptr->comptype == CDI_COMPRESS_ZSTD)
    {
      if (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C
          || streamptr->filetype == CDI_FILETYPE_NCZARR)
        {
          cdfDefVarZstd(streamptr->fileID, ncvarID, streamptr->complevel);
        }
      else
        {
          static bool lwarn = true;
          if (lwarn)
            {
              lwarn = false;
              Warning("SZIP compression is only available for NetCDF4!");
            }
        }
    }
  */
  else if (streamptr->comptype == CDI_COMPRESS_SZIP)
    {
      if (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C
          || streamptr->filetype == CDI_FILETYPE_NCZARR)
        {
          cdfDefVarSzip(streamptr->fileID, ncvarID, xtype2ppb(xtype));
        }
      else
        {
          static bool lwarn = true;
          if (lwarn)
            {
              lwarn = false;
              Warning("SZIP compression is only available for NetCDF4!");
            }
        }
    }
}

static void
cdfDefVarPacking(const stream_t *streamptr, int ncvarID, nc_type xtype, int vlistID, int varID)
{
  //  if ( xtype == NC_BYTE || xtype == NC_SHORT || xtype == NC_INT )
  {
    double addoffset = 0.0, scalefactor = 1.0;
    bool haveAddoffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
    bool haveScalefactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);

    if (haveAddoffset || haveScalefactor)
      {
        nc_type astype = (xtype == NC_FLOAT) ? NC_FLOAT : NC_DOUBLE;
        if ((astype == NC_DOUBLE) && IS_EQUAL(addoffset, (double) ((float) addoffset))
            && IS_EQUAL(scalefactor, (double) ((float) scalefactor)))
          {
            astype = NC_FLOAT;
          }

        cdf_put_att_double(streamptr->fileID, ncvarID, "add_offset", astype, 1, &addoffset);
        cdf_put_att_double(streamptr->fileID, ncvarID, "scale_factor", astype, 1, &scalefactor);
      }
  }
}

static void
cdfAppendCoordinates(int fileID, int ncvarID, char coordinates[CDI_MAX_NAME])
{
  if (ncvarID != CDI_UNDEFID)
    {
      size_t len = strlen(coordinates);
      if (len) coordinates[len++] = ' ';
      cdf_inq_varname(fileID, ncvarID, coordinates + len);
    }
}

static void
cdfDefineCoordinates(const stream_t *streamptr, int ncvarID, int nczvarID, int gridtype, int gridID, int gridindex, int xid,
                     int yid, size_t gridsize, char axis[5], size_t iax)
{
  int fileID = streamptr->fileID;

  if (gridtype != GRID_GENERIC && gridtype != GRID_LONLAT && gridtype != GRID_PROJECTION && gridtype != GRID_CURVILINEAR
      && gridtype != GRID_CHARXY)
    {
      size_t len = strlen(gridNamePtr(gridtype));
      if (len > 0) cdf_put_att_text(fileID, ncvarID, "CDI_grid_type", len, gridNamePtr(gridtype));
    }

  char gmapvarname[CDI_MAX_NAME];
  gmapvarname[0] = 0;
  cdf_get_gmapvarname(gridID, gmapvarname);
  if (gmapvarname[0]) cdf_put_att_text(fileID, ncvarID, "grid_mapping", strlen(gmapvarname), gmapvarname);

  if (gridtype == GRID_GAUSSIAN || gridtype == GRID_GAUSSIAN_REDUCED)
    {
      int numLPE = gridInqNP(gridID);
      if (numLPE > 0) cdf_put_att_int(fileID, ncvarID, "CDI_grid_num_LPE", NC_INT, 1, &numLPE);
    }

  if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      if (ncyvarID != CDI_UNDEFID)
        {
          char name[CDI_MAX_NAME];
          name[0] = 0;
          cdf_inq_varname(fileID, ncyvarID, name);
          size_t len = strlen(name);
          cdf_put_att_text(fileID, ncvarID, "CDI_grid_latitudes", len, name);
        }

      int ncrpvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_RP];
      if (ncrpvarID != CDI_UNDEFID)
        {
          char name[CDI_MAX_NAME];
          name[0] = 0;
          cdf_inq_varname(fileID, ncrpvarID, name);
          size_t len = strlen(name);
          cdf_put_att_text(fileID, ncvarID, "CDI_grid_reduced_points", len, name);
        }
    }

  // define coordinates attribute

  char coordinates[CDI_MAX_NAME];
  coordinates[0] = 0;

  if (nczvarID != CDI_UNDEFID) cdfAppendCoordinates(fileID, nczvarID, coordinates);

  if (gridtype == GRID_TRAJECTORY)
    {
      cdf_put_att_text(fileID, ncvarID, "coordinates", 9, "tlon tlat");
    }
  else if (gridtype == GRID_LONLAT && xid == CDI_UNDEFID && yid == CDI_UNDEFID && gridsize == 1)
    {
      int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      cdfAppendCoordinates(fileID, ncyvarID, coordinates);
      cdfAppendCoordinates(fileID, ncxvarID, coordinates);
    }
  else if (gridtype == GRID_GAUSSIAN_REDUCED)
    {
      /*
      int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      cdfAppendCoordinates(fileID, ncyvarID, coordinates);
      cdfAppendCoordinates(fileID, ncxvarID, coordinates);
      */
    }
  else if (gridtype == GRID_UNSTRUCTURED || gridtype == GRID_CURVILINEAR)
    {
      int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
      int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
      int ncavarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_A];
      // CMOR order: coordinates = "lat lon"
      if (CDI_Coordinates_Lon_Lat)
        {
          cdfAppendCoordinates(fileID, ncxvarID, coordinates);
          cdfAppendCoordinates(fileID, ncyvarID, coordinates);
        }
      else
        {
          cdfAppendCoordinates(fileID, ncyvarID, coordinates);
          cdfAppendCoordinates(fileID, ncxvarID, coordinates);
        }

      if (ncavarID != CDI_UNDEFID)
        {
          char cellarea[CDI_MAX_NAME] = "area: ";
          size_t len = strlen(cellarea);
          cdf_inq_varname(fileID, ncavarID, cellarea + len);
          len = strlen(cellarea);
          cdf_put_att_text(fileID, ncvarID, "cell_measures", len, cellarea);
        }

      if (gridtype == GRID_UNSTRUCTURED)
        {
          int position = gridInqPosition(gridID);
          if (position > 0) cdf_put_att_int(fileID, ncvarID, "number_of_grid_in_reference", NC_INT, 1, &position);
        }
    }
  else if (gridtype == GRID_SPECTRAL || gridtype == GRID_FOURIER)
    {
      axis[iax++] = '-';
      axis[iax++] = '-';
      cdf_put_att_text(fileID, ncvarID, "axis", iax, axis);
      int gridTruncation = gridInqTrunc(gridID);
      cdf_put_att_int(fileID, ncvarID, "truncation", NC_INT, 1, &gridTruncation);
    }
  else if (gridtype == GRID_CHARXY)
    {
      if (gridInqXIsc(gridID))
        {
          int ncxvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_X];
          cdfAppendCoordinates(fileID, ncxvarID, coordinates);
        }
      else if (gridInqYIsc(gridID))
        {
          int ncyvarID = streamptr->ncgrid[gridindex].ncIDs[CDF_VARID_Y];
          cdfAppendCoordinates(fileID, ncyvarID, coordinates);
        }
    }

  size_t len = strlen(coordinates);
  if (len) cdf_put_att_text(fileID, ncvarID, "coordinates", len, coordinates);
}

static size_t
calc_chunksize(size_t chunkSizeLim, size_t size)
{
  static const size_t pageSize = 4096;

  if (size <= chunkSizeLim) return size;

  size_t numChunks = (size / chunkSizeLim) + 1;
  size_t chunkSize = size / numChunks;
  if (chunkSize % pageSize) chunkSize = (chunkSize / pageSize + 1) * pageSize;

  return chunkSize;
}

static const size_t chunkSizeMin = 262144;    // 256k
static const size_t chunkSizeLim = 16777216;  // 16m

size_t
calc_chunksize_y(int chunkType, size_t gridsize, size_t xsize, size_t ysize)
{
  if (chunkType == CDI_CHUNK_AUTO)
    return (gridsize <= chunkSizeMin) ? ysize : chunkSizeMin / xsize;
  else
    return (chunkType == CDI_CHUNK_LINES) ? 1 : ysize;
}

size_t
calc_chunksize_x(int chunkType, int chunkSize, size_t xsize, bool yIsUndefined)
{
  if (chunkType == CDI_CHUNK_AUTO && yIsUndefined)
    return (chunkSize > 0 && chunkSize < (int) xsize) ? (size_t) chunkSize : ((xsize <= chunkSizeMin) ? xsize : chunkSizeMin);
  else
    return calc_chunksize(chunkSizeLim, xsize);
}

static int
cdfDefineDimsAndChunks(const stream_t *streamptr, int varID, int xid, int yid, int zid, size_t gridsize, const int dimorder[3],
                       int dims[4], bool useChunks, size_t chunks[4], char axis[5], size_t *piax)
{
  int fileID = streamptr->fileID;
  int vlistID = streamptr->vlistID;

  size_t iax = *piax;
  int ndims = 0;

  for (int i = 0; i < 4; ++i) chunks[i] = 0;

  size_t xsize = 0, ysize = 0;
  if (xid != CDI_UNDEFID) cdf_inq_dimlen(fileID, xid, &xsize);
  if (yid != CDI_UNDEFID) cdf_inq_dimlen(fileID, yid, &ysize);

  int timetype = vlistInqVarTimetype(vlistID, varID);
  if (vlistHasTime(vlistID) && timetype != TIME_CONSTANT)
    {
      int tid = streamptr->basetime.ncdimid;
      if (tid == CDI_UNDEFID) Error("Internal problem, time undefined!");
      axis[iax++] = 'T';
      chunks[ndims] = 1;
      dims[ndims] = tid;
      ndims++;
    }

  int chunkSize = 0;
  int chunkType = CDI_CHUNK_GRID;
  cdiInqKeyInt(vlistID, varID, CDI_KEY_CHUNKTYPE, &chunkType);
  cdiInqKeyInt(vlistID, varID, CDI_KEY_CHUNKSIZE, &chunkSize);
  if (chunkSize > 0 && yid == CDI_UNDEFID) chunkType = CDI_CHUNK_AUTO;

  if (chunkType == CDI_CHUNK_GRID && gridsize > ChunkSizeLim)
    {
      if (CDI_Debug) fprintf(stderr, "gridsize > %u, changed chunkType to CDI_CHUNK_LINES!\n", ChunkSizeLim);
      chunkType = CDI_CHUNK_LINES;
    }

  for (int id = 0; id < 3; ++id)
    {
      if (dimorder[id] == 3 && zid != CDI_UNDEFID)
        {
          axis[iax++] = 'Z';
          chunks[ndims] = 1;
          dims[ndims] = zid;
          ndims++;
        }
      else if (dimorder[id] == 2 && yid != CDI_UNDEFID)
        {
          chunks[ndims] = calc_chunksize_y(chunkType, gridsize, xsize, ysize);
          dims[ndims] = yid;
          ndims++;
        }
      else if (dimorder[id] == 1 && xid != CDI_UNDEFID)
        {
          chunks[ndims] = calc_chunksize_x(chunkType, chunkSize, xsize, (yid == CDI_UNDEFID));
          dims[ndims] = xid;
          ndims++;
        }
    }

  if (CDI_Debug)
    fprintf(stderr, "useChunks %d chunkType %d chunkSize %d  chunks %zu %zu %zu %zu\n", useChunks, chunkType, chunkSize, chunks[0],
            chunks[1], chunks[2], chunks[3]);

  *piax = iax;
  return ndims;
}

static void
cdfDefineAttrLeveltype(int fileID, int ncvarID, int zaxisID, int zaxistype)
{
  // clang-format off
  if ( zaxistype == ZAXIS_CLOUD_BASE          ||
       zaxistype == ZAXIS_CLOUD_TOP           ||
       zaxistype == ZAXIS_ISOTHERM_ZERO       ||
       zaxistype == ZAXIS_TROPOPAUSE          ||
       zaxistype == ZAXIS_TOA                 ||
       zaxistype == ZAXIS_SEA_BOTTOM          ||
       zaxistype == ZAXIS_LAKE_BOTTOM         ||
       zaxistype == ZAXIS_SEDIMENT_BOTTOM     ||
       zaxistype == ZAXIS_SEDIMENT_BOTTOM_TA  ||
       zaxistype == ZAXIS_SEDIMENT_BOTTOM_TW  ||
       zaxistype == ZAXIS_MIX_LAYER           ||
       zaxistype == ZAXIS_ATMOSPHERE )
    {
      char varname[CDI_MAX_NAME];
      int length = CDI_MAX_NAME;
      cdiInqKeyString(zaxisID, CDI_GLOBAL, CDI_KEY_NAME, varname, &length);
      cdf_put_att_text(fileID, ncvarID, "level_type", strlen(varname), varname);
    }
  // clang-format on
}

static void
cdfDefineAttrEnsemble(int fileID, int ncvarID, int vlistID, int varID)
{
  int perturbationNumber, numberOfForecastsInEnsemble, typeOfEnsembleForecast;
  if (cdiInqKeyInt(vlistID, varID, CDI_KEY_PERTURBATIONNUMBER, &perturbationNumber) == 0)
    cdf_put_att_int(fileID, ncvarID, "realization", NC_INT, 1, &perturbationNumber);
  if (cdiInqKeyInt(vlistID, varID, CDI_KEY_NUMBEROFFORECASTSINENSEMBLE, &numberOfForecastsInEnsemble) == 0)
    cdf_put_att_int(fileID, ncvarID, "ensemble_members", NC_INT, 1, &numberOfForecastsInEnsemble);
  if (cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFENSEMBLEFORECAST, &typeOfEnsembleForecast) == 0)
    cdf_put_att_int(fileID, ncvarID, "forecast_init_type", NC_INT, 1, &typeOfEnsembleForecast);
}

static void
cdfCheckVarname(int fileID, char name[CDI_MAX_NAME])
{
  if (name[0])
    {
      int ncvarID;
      char varname[CDI_MAX_NAME];
      snprintf(varname, sizeof(varname), "%s", name);
      char *varname2 = varname + strlen(varname);
      unsigned iz = 0;

      do
        {
          if (iz) sprintf(varname2, "_%u", iz + 1);

          if (nc_inq_varid(fileID, varname, &ncvarID) != NC_NOERR) break;

          ++iz;
        }
      while (iz <= 99);

      if (iz > 99) Error("Variable name %s already exsist!", name);

      if (!str_is_equal(name, varname))
        Warning("Changed %s entry of variable name '%s' to '%s'!", (iz == 1) ? "double" : "multiple", name, varname);

      strcpy(name, varname);
    }
}

static void
cdfGenVarname(int fileID, char name[CDI_MAX_NAME], int pnum, int pcat, int *pdis, int *pcode)
{
  char varname[CDI_MAX_NAME];

  int code = *pcode;
  if (code < 0) code = -code;
  if (pnum < 0) pnum = -pnum;

  if (*pdis == 255)
    snprintf(varname, sizeof(varname), "var%d", code);
  else
    snprintf(varname, sizeof(varname), "param%d.%d.%d", pnum, pcat, *pdis);

  char *varname2 = varname + strlen(varname);
  int ncvarID;
  unsigned iz = 0;

  do
    {
      if (iz) sprintf(varname2, "_%u", iz + 1);

      if (nc_inq_varid(fileID, varname, &ncvarID) != NC_NOERR) break;

      ++iz;
    }
  while (iz <= 99);

  if (iz > 99) Error("Variable name %s already exsist!", name);

  strcpy(name, varname);
  *pcode = 0;
  *pdis = 255;
}

int
cdfDefVar(stream_t *streamptr, int varID)
{
  if (streamptr->vars[varID].ncvarid != CDI_UNDEFID) return streamptr->vars[varID].ncvarid;

  int fileID = streamptr->fileID;
  if (CDI_Debug) Message("streamID = %d, fileID = %d, varID = %d", streamptr->self, fileID, varID);

  int vlistID = streamptr->vlistID;
  int param = vlistInqVarParam(vlistID, varID);
  int code = vlistInqVarCode(vlistID, varID);
  int pnum, pcat, pdis;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);

  int gridID = vlistInqVarGrid(vlistID, varID);
  SizeType gridsize = gridInqSize(gridID);
  int gridtype = gridInqType(gridID);
  int gridindex = nc_grid_index(streamptr, gridID);
  int xid = (gridtype != GRID_TRAJECTORY) ? streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X] : CDI_UNDEFID;
  int yid = (gridtype != GRID_TRAJECTORY && gridtype != GRID_GAUSSIAN_REDUCED) ? streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y]
                                                                               : CDI_UNDEFID;

  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int zaxistype = zaxisInqType(zaxisID);
  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  int zid = streamptr->zaxisID[zaxisindex];

  int dimorder[3];  // ZYX/321 and ZXY/312
  vlistInqVarDimorder(vlistID, varID, dimorder);
  bool useChunks
      = (gridsize >= 32) && ((dimorder[0] == 3) || (dimorder[1] == 3 && dimorder[2] == 1 && gridsize == gridInqXsize(gridID)));

  if (((dimorder[0] > 0) + (dimorder[1] > 0) + (dimorder[2] > 0))
      < ((xid != CDI_UNDEFID) + (yid != CDI_UNDEFID) + (zid != CDI_UNDEFID)))
    {
      printf("zid=%d  yid=%d  xid=%d\n", zid, yid, xid);
      Error("Internal problem, dimension order missing!");
    }

  size_t iax = 0;
  char axis[5];
  int dims[4];
  size_t chunks[4];
  int ndims = cdfDefineDimsAndChunks(streamptr, varID, xid, yid, zid, gridsize, dimorder, dims, useChunks, chunks, axis, &iax);

  char name[CDI_MAX_NAME];
  int length = CDI_MAX_NAME;
  cdiInqKeyString(vlistID, varID, CDI_KEY_NAME, name, &length);

  char longname[CDI_MAX_NAME];
  vlistInqVarLongname(vlistID, varID, longname);

  char units[CDI_MAX_NAME];
  vlistInqVarUnits(vlistID, varID, units);

  char stdname[CDI_MAX_NAME];
  length = CDI_MAX_NAME;
  cdiInqKeyString(vlistID, varID, CDI_KEY_STDNAME, stdname, &length);

  int tableID = vlistInqVarTable(vlistID, varID);
  if (!name[0]) tableInqEntry(tableID, code, -1, name, longname, units);
  if (name[0])
    cdfCheckVarname(fileID, name);
  else
    cdfGenVarname(fileID, name, pnum, pcat, &pdis, &code);

  int dtype = vlistInqVarDatatype(vlistID, varID);
  const nc_type xtype = cdfDefDatatype(dtype, streamptr);

  if (streamptr->ncmode == 2)
    {
      cdf_redef(fileID);
      streamptr->ncmode = 1;
    }

  int ncvarID = -1;
  cdf_def_var(fileID, name, xtype, ndims, dims, &ncvarID);

#ifdef HAVE_NETCDF4
#ifdef NC_QUANTIZE_BITROUND
  if (xtype == NC_FLOAT || xtype == NC_DOUBLE)
    {
      int nsb = vlistInqVarNSB(vlistID, varID);
      if (nsb > 0) nc_def_var_quantize(fileID, ncvarID, NC_QUANTIZE_BITROUND, nsb);
      // if (nsb > 0) nc_def_var_quantize(fileID, ncvarID, NC_QUANTIZE_GRANULARBR, nsb);
    }
#endif

  if (useChunks
      && (streamptr->filetype == CDI_FILETYPE_NC4 || streamptr->filetype == CDI_FILETYPE_NC4C
          || streamptr->filetype == CDI_FILETYPE_NCZARR))
    cdf_def_var_chunking(fileID, ncvarID, NC_CHUNKED, chunks);
#endif

  if (useChunks) cdfDefVarCompression(streamptr, ncvarID, xtype);
  if (useChunks) cdfDefVarFilter(streamptr, ncvarID);

  if (*stdname) cdf_put_att_text(fileID, ncvarID, "standard_name", strlen(stdname), stdname);
  if (*longname) cdf_put_att_text(fileID, ncvarID, "long_name", strlen(longname), longname);
  if (*units) cdf_put_att_text(fileID, ncvarID, "units", strlen(units), units);

  if (code > 0 && pdis == 255) cdf_put_att_int(fileID, ncvarID, "code", NC_INT, 1, &code);

  if (pdis != 255)
    {
      char paramstr[32];
      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdf_put_att_text(fileID, ncvarID, "param", strlen(paramstr), paramstr);
    }

  if (tableID != CDI_UNDEFID)
    {
      int tablenum = tableInqNum(tableID);
      if (tablenum > 0) cdf_put_att_int(fileID, ncvarID, "table", NC_INT, 1, &tablenum);
    }

  bool zaxisIsScalar = (zid == CDI_UNDEFID) ? (zaxisInqScalar(zaxisID) > 0) : false;
  int nczvarID = (zaxisIsScalar || zaxistype == ZAXIS_CHAR) ? streamptr->nczvarID[zaxisindex] : CDI_UNDEFID;

  cdfDefineCoordinates(streamptr, ncvarID, nczvarID, gridtype, gridID, gridindex, xid, yid, gridsize, axis, iax);

  cdfDefVarPacking(streamptr, ncvarID, xtype, vlistID, varID);

  if (dtype == CDI_DATATYPE_UINT8 && xtype == NC_BYTE)
    {
      int validrange[2] = { 0, 255 };
      cdf_put_att_int(fileID, ncvarID, "valid_range", NC_SHORT, 2, validrange);
      cdf_put_att_text(fileID, ncvarID, "_Unsigned", 4, "true");
    }

  streamptr->vars[varID].ncvarid = ncvarID;

  if (vlistInqVarMissvalUsed(vlistID, varID)) cdfDefVarMissval(streamptr, varID, vlistInqVarDatatype(vlistID, varID), 0);

  if (zid == CDI_UNDEFID) cdfDefineAttrLeveltype(fileID, ncvarID, zaxisID, zaxistype);

  cdfDefineAttrEnsemble(fileID, ncvarID, vlistID, varID);

  // Attribute: cell_methods
  cdfDefineCellMethods(streamptr, vlistID, varID, fileID, ncvarID);

  // Attributes
  cdfDefineAttributes(streamptr->filetype, vlistID, varID, fileID, ncvarID);

  // Institute
  if (vlistInqInstitut(vlistID) == CDI_UNDEFID) cdfDefineInstituteName(vlistID, varID, fileID, ncvarID);

  return ncvarID;
}

void
cdfEndDef(stream_t *streamptr)
{
  cdfDefGlobalAtts(streamptr);

  if (streamptr->accessmode == 0)
    {
      int fileID = streamptr->fileID;
      if (streamptr->ncmode == 2)
        {
          cdf_redef(fileID);
          streamptr->ncmode = 1;
        }

      int nvars = streamptr->nvars;
      for (int varID = 0; varID < nvars; ++varID) cdfDefVar(streamptr, varID);

      if (streamptr->ncmode != 2)
        {
          if (CDI_Netcdf_Hdr_Pad == 0UL)
            cdf_enddef(fileID, streamptr->self);
          else
            cdf__enddef(fileID, streamptr->self, CDI_Netcdf_Hdr_Pad);
          streamptr->ncmode = 2;
        }

      streamptr->accessmode = 1;
    }
}

static void
cdfWriteGridTraj(stream_t *streamptr, int gridID)
{
  int gridindex = nc_grid_index(streamptr, gridID);
  int lonID = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
  int latID = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
  size_t index = (size_t) streamptr->curTsID;

  double xlon = gridInqXval(gridID, 0);
  double xlat = gridInqYval(gridID, 0);

  cdf_put_var1_double(streamptr->fileID, lonID, &index, &xlon);
  cdf_put_var1_double(streamptr->fileID, latID, &index, &xlat);
}

static void
cdf_write_var_data(int fileID, int vlistID, int varID, int ncvarID, int dtype, size_t nvals, size_t xsize, size_t ysize,
                   bool swapxy, size_t *start, size_t *count, int memtype, const void *data, size_t nmiss)
{
  const double *pdata_dp = (const double *) data;
  double *mdata_dp = NULL;
  double *sdata_dp = NULL;
  const float *pdata_sp = (const float *) data;
  float *mdata_sp = NULL;
  float *sdata_sp = NULL;

  /*  if ( dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_INT16 || dtype == CDI_DATATYPE_INT32 ) */
  {
    double missval = vlistInqVarMissval(vlistID, varID);
    double addoffset = 0.0, scalefactor = 1.0;
    bool haveAddoffset = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_ADDOFFSET, &addoffset) == CDI_NOERR);
    bool haveScalefactor = (cdiInqKeyFloat(vlistID, varID, CDI_KEY_SCALEFACTOR, &scalefactor) == CDI_NOERR);

    if (haveAddoffset || haveScalefactor)
      {
        if (memtype == MEMTYPE_FLOAT)
          {
            mdata_sp = (float *) Malloc(nvals * sizeof(float));
            memcpy(mdata_sp, pdata_sp, nvals * sizeof(float));
            pdata_sp = mdata_sp;

            if (nmiss > 0)
              {
                for (size_t i = 0; i < nvals; ++i)
                  {
                    double temp = mdata_sp[i];
                    if (!DBL_IS_EQUAL(temp, (float) missval))
                      {
                        if (haveAddoffset) temp -= addoffset;
                        if (haveScalefactor) temp /= scalefactor;
                        mdata_sp[i] = (float) temp;
                      }
                  }
              }
            else
              {
                for (size_t i = 0; i < nvals; ++i)
                  {
                    double temp = mdata_sp[i];
                    if (haveAddoffset) temp -= addoffset;
                    if (haveScalefactor) temp /= scalefactor;
                    mdata_sp[i] = (float) temp;
                  }
              }
          }
        else
          {
            mdata_dp = (double *) Malloc(nvals * sizeof(double));
            memcpy(mdata_dp, pdata_dp, nvals * sizeof(double));
            pdata_dp = mdata_dp;

            if (nmiss > 0)
              {
                for (size_t i = 0; i < nvals; ++i)
                  {
                    if (!DBL_IS_EQUAL(mdata_dp[i], missval))
                      {
                        if (haveAddoffset) mdata_dp[i] -= addoffset;
                        if (haveScalefactor) mdata_dp[i] /= scalefactor;
                      }
                  }
              }
            else
              {
                for (size_t i = 0; i < nvals; ++i)
                  {
                    if (haveAddoffset) mdata_dp[i] -= addoffset;
                    if (haveScalefactor) mdata_dp[i] /= scalefactor;
                  }
              }
          }
      }

    if (dtype == CDI_DATATYPE_UINT8 || dtype == CDI_DATATYPE_INT8 || dtype == CDI_DATATYPE_UINT16 || dtype == CDI_DATATYPE_INT16
        || dtype == CDI_DATATYPE_UINT32 || dtype == CDI_DATATYPE_INT32)
      {
        if (memtype == MEMTYPE_FLOAT)
          {
            if (mdata_sp == NULL)
              {
                mdata_sp = (float *) Malloc(nvals * sizeof(float));
                memcpy(mdata_sp, pdata_sp, nvals * sizeof(float));
                pdata_sp = mdata_sp;
              }

            for (size_t i = 0; i < nvals; ++i) mdata_sp[i] = roundf(mdata_sp[i]);

            if (dtype == CDI_DATATYPE_UINT8)
              {
                nc_type xtype;
                cdf_inq_vartype(fileID, ncvarID, &xtype);
                if (xtype == NC_BYTE)
                  {
                    for (size_t i = 0; i < nvals; ++i)
                      if (mdata_sp[i] > 127) mdata_sp[i] -= 256;
                  }
              }
          }
        else
          {
            if (mdata_dp == NULL)
              {
                mdata_dp = (double *) Malloc(nvals * sizeof(double));
                memcpy(mdata_dp, pdata_dp, nvals * sizeof(double));
                pdata_dp = mdata_dp;
              }

            for (size_t i = 0; i < nvals; ++i) mdata_dp[i] = round(mdata_dp[i]);

            if (dtype == CDI_DATATYPE_UINT8)
              {
                nc_type xtype;
                cdf_inq_vartype(fileID, ncvarID, &xtype);
                if (xtype == NC_BYTE)
                  {
                    for (size_t i = 0; i < nvals; ++i)
                      if (mdata_dp[i] > 127) mdata_dp[i] -= 256;
                  }
              }
          }
      }

    if (CDF_Debug)
      {
        double fmin = 1.0e200;
        double fmax = -1.0e200;
        if (memtype == MEMTYPE_FLOAT)
          {
            for (size_t i = 0; i < nvals; ++i)
              {
                if (!DBL_IS_EQUAL(pdata_sp[i], (float) missval))
                  {
                    if (pdata_sp[i] < fmin) fmin = pdata_sp[i];
                    if (pdata_sp[i] > fmax) fmax = pdata_sp[i];
                  }
              }
          }
        else
          {
            for (size_t i = 0; i < nvals; ++i)
              {
                if (!DBL_IS_EQUAL(pdata_dp[i], missval))
                  {
                    if (pdata_dp[i] < fmin) fmin = pdata_dp[i];
                    if (pdata_dp[i] > fmax) fmax = pdata_dp[i];
                  }
              }
          }

        Message("nvals = %zu, nmiss = %d, missval = %g, minval = %g, maxval = %g", nvals, nmiss, missval, fmin, fmax);
      }
  }

  if (swapxy)  // implemented only for cdf_write_var_slice()
    {
      size_t gridsize = xsize * ysize;
      if (memtype == MEMTYPE_FLOAT)
        {
          sdata_sp = (float *) Malloc(gridsize * sizeof(float));
          for (size_t j = 0; j < ysize; ++j)
            for (size_t i = 0; i < xsize; ++i) sdata_sp[i * ysize + j] = pdata_sp[j * xsize + i];
          pdata_sp = sdata_sp;
        }
      else
        {
          sdata_dp = (double *) Malloc(gridsize * sizeof(double));
          for (size_t j = 0; j < ysize; ++j)
            for (size_t i = 0; i < xsize; ++i) sdata_dp[i * ysize + j] = pdata_dp[j * xsize + i];
          pdata_dp = sdata_dp;
        }
    }

  if (dtype == CDI_DATATYPE_CPX32 || dtype == CDI_DATATYPE_CPX64)
    {
      void *pdata = (memtype == MEMTYPE_FLOAT) ? (void *) pdata_sp : (void *) pdata_dp;
      float *cdata_sp = NULL;
      double *cdata_dp = NULL;
      if (memtype == MEMTYPE_FLOAT && dtype == CDI_DATATYPE_CPX64)
        {
          cdata_dp = (double *) Malloc(2 * nvals * sizeof(double));
          for (size_t i = 0; i < nvals; ++i)
            {
              cdata_dp[2 * i] = (double) (pdata_sp[2 * i]);
              cdata_dp[2 * i + 1] = (double) (pdata_sp[2 * i + 1]);
            }
          pdata = cdata_dp;
        }
      else if (memtype == MEMTYPE_DOUBLE && dtype == CDI_DATATYPE_CPX32)
        {
          cdata_sp = (float *) Malloc(2 * nvals * sizeof(float));
          for (size_t i = 0; i < nvals; ++i)
            {
              cdata_sp[2 * i] = (float) (pdata_dp[2 * i]);
              cdata_sp[2 * i + 1] = (float) (pdata_dp[2 * i + 1]);
            }
          pdata = cdata_sp;
        }

      cdf_put_vara(fileID, ncvarID, start, count, pdata);
      if (cdata_sp) Free(cdata_sp);
      if (cdata_dp) Free(cdata_dp);
    }
  else
    {
      if (memtype == MEMTYPE_FLOAT)
        cdf_put_vara_float(fileID, ncvarID, start, count, pdata_sp);
      else
        cdf_put_vara_double(fileID, ncvarID, start, count, pdata_dp);
    }

  if (mdata_dp) Free(mdata_dp);
  if (sdata_dp) Free(sdata_dp);
  if (mdata_sp) Free(mdata_sp);
  if (sdata_sp) Free(sdata_sp);
}

static void
cdfGetXYZid(stream_t *streamptr, int gridID, int zaxisID, int *xid, int *yid, int *zid)
{
  *xid = CDI_UNDEFID;
  *yid = CDI_UNDEFID;

  int gridtype = gridInqType(gridID);
  if (gridtype == GRID_TRAJECTORY)
    {
      cdfWriteGridTraj(streamptr, gridID);
    }
  else
    {
      int gridindex = nc_grid_index(streamptr, gridID);
      *xid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_X];
      if (gridtype != GRID_GAUSSIAN_REDUCED) *yid = streamptr->ncgrid[gridindex].ncIDs[CDF_DIMID_Y];
    }

  int vlistID = streamptr->vlistID;
  int zaxisindex = vlistZaxisIndex(vlistID, zaxisID);
  *zid = streamptr->zaxisID[zaxisindex];
}

static void
cdfDefineStartAndCount(stream_t *streamptr, int varID, int xid, int yid, int zid, size_t start[5], size_t count[5], size_t *xsize,
                       size_t *ysize)
{
  size_t ndims = 0;
  *xsize = 0;
  *ysize = 0;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  const long ntsteps = streamptr->ntsteps;
  if (CDI_Debug) Message("ntsteps = %ld", ntsteps);

  int timetype = vlistInqVarTimetype(vlistID, varID);

  if (vlistHasTime(vlistID) && timetype != TIME_CONSTANT)
    {
      start[ndims] = (size_t) ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  if (zid != CDI_UNDEFID)
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      start[ndims] = 0;
      count[ndims] = (size_t) zaxisInqSize(zaxisID);
      ndims++;
    }

  if (yid != CDI_UNDEFID)
    {
      start[ndims] = 0;
      size_t size;
      cdf_inq_dimlen(fileID, yid, &size);
      /*      count[ndims] = gridInqYsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if (xid != CDI_UNDEFID)
    {
      start[ndims] = 0;
      size_t size;
      cdf_inq_dimlen(fileID, xid, &size);
      /*      count[ndims] = gridInqXsize(gridID); */
      count[ndims] = size;
      ndims++;
    }

  if (CDI_Debug)
    for (size_t idim = 0; idim < ndims; ++idim) Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);
}

void
cdf_write_var(stream_t *streamptr, int varID, int memtype, const void *data, size_t nmiss)
{
  if (streamptr->accessmode == 0) cdfEndDef(streamptr);

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  int ncvarID = cdfDefVar(streamptr, varID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  int xid, yid, zid;
  cdfGetXYZid(streamptr, gridID, zaxisID, &xid, &yid, &zid);

  size_t xsize, ysize;
  size_t start[5], count[5];
  cdfDefineStartAndCount(streamptr, varID, xid, yid, zid, start, count, &xsize, &ysize);

  if (streamptr->ncmode == 1)
    {
      cdf_enddef(fileID, streamptr->self);
      streamptr->ncmode = 2;
    }

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if (nmiss > 0) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = gridInqSize(gridID) * (size_t) (zaxisInqSize(zaxisID));

  bool swapxy = false;
  cdf_write_var_data(fileID, vlistID, varID, ncvarID, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}

static void
cdfDefineStartAndCountChunk(stream_t *streamptr, const int rect[][2], int varID, int xid, int yid, int zid, size_t start[5],
                            size_t count[5], size_t *xsize, size_t *ysize)
{
  size_t ndims = 0;
  *xsize = 0;
  *ysize = 0;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  const long ntsteps = streamptr->ntsteps;
  if (CDI_Debug) Message("ntsteps = %ld", ntsteps);

  int timetype = vlistInqVarTimetype(vlistID, varID);

  if (vlistHasTime(vlistID) && timetype != TIME_CONSTANT)
    {
      start[ndims] = (size_t) ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  if (zid != CDI_UNDEFID)
    {
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int size = zaxisInqSize(zaxisID);
      xassert(rect[2][0] >= 0 && rect[2][0] <= rect[2][1] + 1 && rect[2][1] <= size);
      start[ndims] = (size_t) rect[2][0];
      count[ndims] = rect[2][1] < 0 ? (size_t) 0 : (size_t) rect[2][1] - (size_t) rect[2][0] + 1;
      ndims++;
    }

  if (yid != CDI_UNDEFID)
    {
      size_t size;
      cdf_inq_dimlen(fileID, yid, &size);
      xassert(rect[1][0] >= 0 && rect[1][0] <= rect[1][1] + 1 && (rect[1][1] < 0 || (size_t) rect[1][1] <= size));
      start[ndims] = (size_t) rect[1][0];
      count[ndims] = rect[1][1] < 0 ? (size_t) 0 : ((size_t) rect[1][1] - (size_t) rect[1][0] + 1);
      ndims++;
    }

  if (xid != CDI_UNDEFID)
    {
      size_t size;
      cdf_inq_dimlen(fileID, xid, &size);
      xassert(rect[0][0] >= 0 && rect[0][0] <= rect[0][1] + 1 && (rect[0][1] < 0 || (size_t) rect[0][1] <= size));
      start[ndims] = (size_t) rect[0][0];
      count[ndims] = rect[0][1] < 0 ? (size_t) 0 : (size_t) rect[0][1] - (size_t) rect[0][0] + 1;
      ndims++;
    }

  if (CDI_Debug)
    for (size_t idim = 0; idim < ndims; ++idim) Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);
}

void
cdf_write_var_chunk(stream_t *streamptr, int varID, int memtype, const int rect[][2], const void *data, size_t nmiss)
{
  if (streamptr->accessmode == 0) cdfEndDef(streamptr);

  int streamID = streamptr->self;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  int vlistID = streamInqVlist(streamID);
  int fileID = streamInqFileID(streamID);

  int ncvarID = cdfDefVar(streamptr, varID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  int xid, yid, zid;
  cdfGetXYZid(streamptr, gridID, zaxisID, &xid, &yid, &zid);

  size_t xsize, ysize;
  size_t start[5], count[5];
  cdfDefineStartAndCountChunk(streamptr, rect, varID, xid, yid, zid, start, count, &xsize, &ysize);

  if (streamptr->ncmode == 1)
    {
      cdf_enddef(fileID, streamptr->self);
      streamptr->ncmode = 2;
    }

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if (nmiss > 0) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = gridInqSize(gridID) * (size_t) (zaxisInqSize(zaxisID));

  bool swapxy = false;
  cdf_write_var_data(fileID, vlistID, varID, ncvarID, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}

static void
cdfDefineStartAndCountSlice(stream_t *streamptr, int varID, int levelID, int dimorder[3], int xid, int yid, int zid,
                            size_t start[5], size_t count[5], size_t *xsize, size_t *ysize)
{
  size_t ndims = 0;
  *xsize = 0;
  *ysize = 0;

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  const long ntsteps = streamptr->ntsteps;
  if (CDI_Debug) Message("ntsteps = %ld", ntsteps);

  int timetype = vlistInqVarTimetype(vlistID, varID);

  if (vlistHasTime(vlistID) && timetype != TIME_CONSTANT)
    {
      start[ndims] = (size_t) ntsteps - 1;
      count[ndims] = 1;
      ndims++;
    }

  for (int id = 0; id < 3; ++id)
    {
      if (dimorder[id] == 3 && zid != CDI_UNDEFID)
        {
          start[ndims] = (size_t) levelID;
          count[ndims] = 1;
          ndims++;
        }
      else if (dimorder[id] == 2 && yid != CDI_UNDEFID)
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, yid, ysize);
          count[ndims] = *ysize;
          ndims++;
        }
      else if (dimorder[id] == 1 && xid != CDI_UNDEFID)
        {
          start[ndims] = 0;
          cdf_inq_dimlen(fileID, xid, xsize);
          count[ndims] = *xsize;
          ndims++;
        }
    }

  if (CDI_Debug)
    for (size_t idim = 0; idim < ndims; ++idim) Message("dim = %d  start = %d  count = %d", idim, start[idim], count[idim]);
}

void
cdf_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, size_t nmiss)
{
  if (streamptr->accessmode == 0) cdfEndDef(streamptr);

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamptr->self, varID);

  int vlistID = streamptr->vlistID;
  int fileID = streamptr->fileID;

  int ncvarID = cdfDefVar(streamptr, varID);

  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);

  int xid, yid, zid;
  cdfGetXYZid(streamptr, gridID, zaxisID, &xid, &yid, &zid);

  int dimorder[3];
  vlistInqVarDimorder(vlistID, varID, dimorder);
  bool swapxy = (dimorder[2] == 2 || dimorder[0] == 1) && xid != CDI_UNDEFID && yid != CDI_UNDEFID;

  size_t xsize, ysize;
  size_t start[5], count[5];
  cdfDefineStartAndCountSlice(streamptr, varID, levelID, dimorder, xid, yid, zid, start, count, &xsize, &ysize);

  int dtype = vlistInqVarDatatype(vlistID, varID);

  if (nmiss > 0) cdfDefVarMissval(streamptr, varID, dtype, 1);

  size_t nvals = gridInqSize(gridID);

  cdf_write_var_data(fileID, vlistID, varID, ncvarID, dtype, nvals, xsize, ysize, swapxy, start, count, memtype, data, nmiss);
}

void
cdf_write_record(stream_t *streamptr, int memtype, const void *data, size_t nmiss)
{
  int varID = streamptr->record->varID;
  int levelID = streamptr->record->levelID;
  cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
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
