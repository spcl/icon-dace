#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBGRIB

#ifdef HAVE_LIBFDB5
#include "cdi_fdb.h"
#endif

#ifdef HAVE_ACROSS
#include "cdi_across.h"
#endif

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "gribapi.h"
#include "stream_cgribex.h"
#include "stream_grb.h"
#include "stream_gribapi.h"
#include "file.h"
#include "cgribex.h" /* gribZip gribGetZip gribGinfo */
#include "namespace.h"

static size_t
grb_encode(int filetype, int memType, int varID, int levelID, int vlistID, int gridID, int zaxisID, CdiDateTime vDateTime,
           int tsteptype, int numavg, size_t datasize, const void *data, size_t nmiss, void **gribbuffer, int comptype,
           void *gribContainers)
{
  size_t nbytes = 0;

#ifdef HAVE_LIBCGRIBEX
  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1)
    {
      size_t gribbuffersize = datasize * 4 + 3000;
      *gribbuffer = Malloc(gribbuffersize);

      nbytes = cgribexEncode(memType, varID, levelID, vlistID, gridID, zaxisID, vDateTime, tsteptype, numavg, datasize, data, nmiss,
                             *gribbuffer, gribbuffersize);
    }
  else
#endif
#ifdef HAVE_LIBGRIB_API
    {
#ifdef GRIBCONTAINER2D
      void *gribContainer = (void *) &((gribContainer_t **) gribContainers)[varID][levelID];
#else
      void *gribContainer = (void *) &((gribContainer_t *) gribContainers)[varID];
#endif
      const void *datap = data;
      int memTypeX = /*have_gribapi_float_interface() ? memType :*/ MEMTYPE_DOUBLE;
      if (/*!have_gribapi_float_interface() &&*/ memType == MEMTYPE_FLOAT)
        {
          const float *dataf = (const float *) data;
          double *datad = (double *) Malloc(datasize * sizeof(double));
          for (size_t i = 0; i < datasize; ++i) datad[i] = (double) dataf[i];
          datap = (const void *) datad;
        }

      size_t gribbuffersize;
      nbytes = gribapiEncode(memTypeX, varID, levelID, vlistID, gridID, zaxisID, vDateTime, tsteptype, numavg, (long) datasize,
                             datap, nmiss, gribbuffer, &gribbuffersize, comptype, gribContainer);

      if (/*!have_gribapi_float_interface() &&*/ memType == MEMTYPE_FLOAT) Free((void *) datap);
    }
#else
  {
    Error("ecCodes support not compiled in!");
    (void) gribContainers;
    (void) comptype;
  }
#endif

  return nbytes;
}

static size_t
grbSzip(int filetype, void *gribbuffer, size_t gribbuffersize)
{
  size_t buffersize = gribbuffersize + 1000;  // compressed record can be greater than source record
  void *buffer = Malloc(buffersize);

  // memcpy(buffer, gribbuffer, gribbuffersize);

  size_t nbytes = 0;
  if (filetype == CDI_FILETYPE_GRB)
    {
      nbytes = (size_t) gribZip((unsigned char *) gribbuffer, (long) gribbuffersize, (unsigned char *) buffer, (long) buffersize);
    }
  else
    {
      static int lszip_warn = 1;
      if (lszip_warn) Warning("Szip compression of GRIB2 records not implemented!");
      lszip_warn = 0;
      nbytes = gribbuffersize;
    }

  Free(buffer);

  return nbytes;
}

typedef struct
{
  char date[16];
  char time[8];
  char param[8];
  char levtype[8];
  char levelist[8];
} FDB_Keys;

void
cdi_fdb_store(void *fdbHandle, const char *filename, void *gribbuffer, size_t nbytes, const FDB_Keys *fdbKeys)
{
#ifdef HAVE_LIBFDB5
  size_t len = strlen(filename);
  if (len == 4) Error("FDB keys missing!");

  KeyValueEntry keyValue;
  keyValue.item = NULL;
  decode_fdbitem(filename + 4, &keyValue);

  if (keyValue.numKeys == 0) Error("FDB keys missing!");

  const char *class = NULL;
  const char *stream = NULL;
  const char *expver = NULL;
  for (int i = 0; i < keyValue.numKeys; i++)
    {
      // clang-format off
      if      (str_is_equal(keyValue.keys[i], "class"))  class = keyValue.values[i];
      else if (str_is_equal(keyValue.keys[i], "stream")) stream = keyValue.values[i];
      else if (str_is_equal(keyValue.keys[i], "expver")) expver = keyValue.values[i];
      else Error("Unsupported FDB parameter: %s=%s", keyValue.keys[i], keyValue.values[i]);
      // clang-format on
    }

  if (!class) Error("FDB parameter <class> undefined!");
  if (!stream) Error("FDB parameter <stream> undefined!");
  if (!expver) Error("FDB parameter <expver> undefined!");

  if (CDI_Debug)
    {
      printf("{class=%s,expver=%s,stream=%s,date=%s,time=%s,domain=g}", class, expver, stream, fdbKeys->date, fdbKeys->time);
      printf("{type=an,levtype=%s}{step=0,", fdbKeys->levtype);
      if (fdbKeys->levelist[0]) printf("levelist=%s,", fdbKeys->levelist);
      printf("param=%s},length=%zu\n", fdbKeys->param, nbytes);
    }

  fdb_key_t *key;
  fdb_new_key(&key);
  fdb_key_add(key, "class", class);
  fdb_key_add(key, "expver", expver);
  fdb_key_add(key, "stream", stream);
  fdb_key_add(key, "domain", "g");
  fdb_key_add(key, "date", fdbKeys->date);
  fdb_key_add(key, "time", fdbKeys->time);
  fdb_key_add(key, "type", "an");
  fdb_key_add(key, "levtype", fdbKeys->levtype);
  fdb_key_add(key, "step", "0");
  fdb_key_add(key, "param", fdbKeys->param);
  if (fdbKeys->levelist[0]) fdb_key_add(key, "levelist", fdbKeys->levelist);

  fdb_archive(fdbHandle, key, gribbuffer, nbytes);
  // alternative: write to tmpfile and use C++ code from fdb_write

  fdb_delete_key(key);
#endif
}

static int
map_grib2_param(int pnum, int pcat, int pdis)
{
  // clang-format off
  if      (pnum ==  8 && pcat == 2 && pdis == 0) return 135;
  else if (pnum ==  0 && pcat == 0 && pdis == 0) return 130;
  else if (pnum ==  0 && pcat == 1 && pdis == 0) return 133;
  else if (pnum == 83 && pcat == 1 && pdis == 0) return 246;
  else if (pnum == 84 && pcat == 1 && pdis == 0) return 247;
  else if (pnum == 85 && pcat == 1 && pdis == 0) return  75;
  else if (pnum == 86 && pcat == 1 && pdis == 0) return  76;
  else if (pnum ==  2 && pcat == 2 && pdis == 0) return 131;
  else if (pnum ==  3 && pcat == 2 && pdis == 0) return 132;
  else if (pnum == 25 && pcat == 3 && pdis == 0) return 152;
  else if (pnum ==  4 && pcat == 3 && pdis == 0) return 129;
  // clang-format on

  return -1;
}

static int
get_fdb_param(int cdiParam)
{
  int pnum, pcat, pdis;
  cdiDecodeParam(cdiParam, &pnum, &pcat, &pdis);
  if (pnum < 0) pnum = -pnum;
  if (pnum > 255) pnum = pnum % 256;

  int fdbParam = (pdis == 255) ? pnum : map_grib2_param(pnum, pcat, pdis);

  return fdbParam;
}

static void
fillup_gribbuffer(size_t nbytes, unsigned char *gribbuffer)
{
  while (nbytes & 7) gribbuffer[nbytes++] = 0;
}

void
grbCopyRecord(stream_t *streamptr2, stream_t *streamptr1)
{
  int filetype = streamptr1->filetype;
  int fileID1 = streamptr1->fileID;
  int fileID2 = streamptr2->fileID;
  int tsID = streamptr1->curTsID;
  int vrecID = streamptr1->tsteps[tsID].curRecID;
  int recID = streamptr1->tsteps[tsID].recIDs[vrecID];
  const record_t *record = &streamptr1->tsteps[tsID].records[recID];
  off_t recpos = record->position;
  size_t recsize = record->size;

  void *gribbuffer = NULL;

  if (streamptr1->protocol == CDI_PROTOCOL_FDB)
    {
#ifdef HAVE_LIBFDB5
      void *fdbItem = streamptr1->tsteps[tsID].records[recID].fdbItem;
      if (!fdbItem) Error("fdbItem not available!");

      size_t buffersize = 0;
      recsize = fdb_read_record(streamptr1->protocolData, fdbItem, &buffersize, &gribbuffer);

      // round up recsize to next multiple of 8
      size_t gribbuffersize = ((recsize + 7U) & ~7U);

      gribbuffer = (unsigned char *) Realloc(gribbuffer, gribbuffersize);
#endif
    }
  else
    {
      fileSetPos(fileID1, recpos, SEEK_SET);

      // round up recsize to next multiple of 8
      size_t gribbuffersize = ((recsize + 7U) & ~7U);

      gribbuffer = (unsigned char *) Malloc(gribbuffersize);

      if (fileRead(fileID1, gribbuffer, recsize) != recsize) Error("Could not read GRIB record for copying!");
    }

  size_t nbytes = recsize;

#ifdef HAVE_LIBCGRIBEX
  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1)
    {
      if (cdiGribChangeParameterID.active)
        {
          // Even if you are stream-copy records you might need to change a bit of grib-header !
          void *gh = cgribex_handle_new_from_meassage((void *) gribbuffer, recsize);
          cgribexChangeParameterIdentification(gh, cdiGribChangeParameterID.code, cdiGribChangeParameterID.ltype,
                                               cdiGribChangeParameterID.lev);
          cgribex_handle_delete(gh);
          cdiGribChangeParameterID.active = false;  // after grib attributes have been changed turn it off again
        }
    }
  else
#endif
#ifdef HAVE_LIBGRIB_API
    {
      if (cdiGribChangeParameterID.active)
        {
          // Even if you are stream-copy records you might need to change a bit of grib-header !
          grib_handle *gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);
          gribapiChangeParameterIdentification(gh, cdiGribChangeParameterID.code, cdiGribChangeParameterID.ltype,
                                               cdiGribChangeParameterID.lev);
          grib_handle_delete(gh);
          cdiGribChangeParameterID.active = false;  // after grib attributes have been changed turn it off again
        }
    }
#else
  {
    Error("ecCodes support not compiled in!");
  }
#endif

#ifdef HAVE_LIBGRIB_API
#ifdef HIRLAM_EXTENSIONS
  // Even if you are stream-copy records you might need to change a bit of grib-header !

  if (cdiGribDataScanningMode.active)
    // allowed modes: <0, 64, 96>; Default is 64
    // This will overrule the old scanning mode of the given grid
    {
      grib_handle *gh = grib_handle_new_from_message(NULL, (void *) gribbuffer, recsize);

      int scanModeIN = gribapiGetScanningMode(gh);

      grib_handle_delete(gh);

      if (cdiDebugExt >= 20)
        Message("Change GribDataScanningMode => %d (scanModeIN = %d)", cdiGribDataScanningMode.value, scanModeIN);

      if (scanModeIN != cdiGribDataScanningMode.value)
        {
          size_t nmiss = 0;

          int vlistID = streamptr1->vlistID;
          int varID = record->varID;
          int levelID = record->levelID;
          int gridID = vlistInqVarGrid(vlistID, varID);

          size_t gridsize = gridInqSize(gridID);
          if (vlistNumber(vlistID) != CDI_REAL) gridsize *= 2;
          double *data = (double *) Malloc(gridsize * sizeof(double));

          if (cdiDebugExt >= 20) Message(" processing varID %d; levelID %d", varID, levelID);

          grb_write_var_slice(streamptr2, varID, levelID, MEMTYPE_DOUBLE, (const void *) data, nmiss);

          free(data);
        }
    }
#endif  // HIRLAM_EXTENSIONS
#endif

  if (filetype == CDI_FILETYPE_GRB)
    {
      size_t unzipsize;
      int izip = gribGetZip(recsize, (unsigned char *) gribbuffer, &unzipsize);

      if (izip == 0 && (streamptr2->comptype == CDI_COMPRESS_SZIP || streamptr2->comptype == CDI_COMPRESS_AEC))
        nbytes = grbSzip(filetype, gribbuffer, nbytes);
    }

  fillup_gribbuffer(nbytes, (unsigned char *) gribbuffer);

  if (streamptr2->protocol == CDI_PROTOCOL_FDB)
    {
      int vlistID = streamptr1->vlistID;
      int varID = record->varID;
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int zaxisType = zaxisInqType(zaxisID);
      CdiDateTime vDateTime = streamptr1->tsteps[tsID].taxis.vDateTime;

      FDB_Keys fdbKeys;
      snprintf(fdbKeys.date, sizeof(fdbKeys.date), "%d", (int) cdiDate_get(vDateTime.date));
      snprintf(fdbKeys.time, sizeof(fdbKeys.time), "%04d", (short) (cdiTime_get(vDateTime.time) / 100));
      snprintf(fdbKeys.param, sizeof(fdbKeys.param), "%d", get_fdb_param(record->param));
      bool isML = (zaxisType == ZAXIS_HYBRID || zaxisType == ZAXIS_HYBRID_HALF);
      snprintf(fdbKeys.levtype, sizeof(fdbKeys.levtype), "%s", isML ? "ml" : "sfc");
      fdbKeys.levelist[0] = 0;
      if (isML) snprintf(fdbKeys.levelist, sizeof(fdbKeys.levelist), "%d", isML ? record->ilevel : 0);

#ifdef HAVE_LIBFDB5
      cdi_fdb_store(streamptr2->protocolData, streamptr2->filename, gribbuffer, nbytes, &fdbKeys);
#endif
    }
  else
    {
      size_t nwrite = fileWrite(fileID2, gribbuffer, nbytes);
      if (nwrite != nbytes) SysError("Could not write record for copying!");
    }

  Free(gribbuffer);
}

void
grb_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, size_t nmiss)
{
  int filetype = streamptr->filetype;
  int fileID = streamptr->fileID;
  int vlistID = streamptr->vlistID;
  int gridID = vlistInqVarGrid(vlistID, varID);
  int zaxisID = vlistInqVarZaxis(vlistID, varID);
  int tsteptype = vlistInqVarTsteptype(vlistID, varID);
  int tsID = streamptr->curTsID;
  CdiDateTime vDateTime = streamptr->tsteps[tsID].taxis.vDateTime;
  int numavg = (tsteptype == TSTEP_AVG) ? streamptr->tsteps[tsID].taxis.numavg : 0;
  int comptype = streamptr->comptype;

  if (CDI_Debug) Message("gridID = %d zaxisID = %d", gridID, zaxisID);

  size_t datasize = gridInqSize(gridID);

  if (comptype != CDI_COMPRESS_JPEG && comptype != CDI_COMPRESS_SZIP && comptype != CDI_COMPRESS_AEC) comptype = CDI_COMPRESS_NONE;

  if (filetype == CDI_FILETYPE_GRB && !CDI_gribapi_grib1 && comptype == CDI_COMPRESS_JPEG)
    {
      static bool ljpeg_warn = true;
      if (ljpeg_warn) Warning("JPEG compression of GRIB1 records not available!");
      ljpeg_warn = false;
    }

  void *gribbuffer = NULL;
  size_t nbytes = grb_encode(filetype, memtype, varID, levelID, vlistID, gridID, zaxisID, vDateTime, tsteptype, numavg, datasize,
                             data, nmiss, &gribbuffer, comptype, streamptr->gribContainers);

  if (filetype == CDI_FILETYPE_GRB && (comptype == CDI_COMPRESS_SZIP || comptype == CDI_COMPRESS_AEC))
    nbytes = grbSzip(filetype, gribbuffer, nbytes);

  switch (streamptr->protocol)
    {
    case CDI_PROTOCOL_ACROSS:
      {
#ifdef HAVE_ACROSS
        if (across_write_grib_message(streamptr, gribbuffer, nbytes)) SysError("Failed to write GRIB slice!");
#endif
      }
      break;

    case CDI_PROTOCOL_FDB:
      {
#ifdef HAVE_LIBFDB5
        int zaxisType = zaxisInqType(zaxisID);
        double level = zaxisInqLevel(zaxisID, levelID);

        FDB_Keys fdbKeys;
        snprintf(fdbKeys.date, sizeof(fdbKeys.date), "%d", (int) cdiDate_get(vDateTime.date));
        snprintf(fdbKeys.time, sizeof(fdbKeys.time), "%04d", (short) (cdiTime_get(vDateTime.time) / 100));
        snprintf(fdbKeys.param, sizeof(fdbKeys.param), "%d", get_fdb_param(vlistInqVarParam(vlistID, varID)));
        bool isML = (zaxisType == ZAXIS_HYBRID || zaxisType == ZAXIS_HYBRID_HALF);
        snprintf(fdbKeys.levtype, sizeof(fdbKeys.levtype), "%s", isML ? "ml" : "sfc");
        fdbKeys.levelist[0] = 0;
        int ilevel = (isML) ? (int) level : 0;
        if (isML) snprintf(fdbKeys.levelist, sizeof(fdbKeys.levelist), "%d", isML ? ilevel : 0);

        cdi_fdb_store(streamptr->protocolData, streamptr->filename, gribbuffer, nbytes, &fdbKeys);
#endif
      }
      break;

    case CDI_PROTOCOL_OTHER:
    case CDI_PROTOCOL_FILE:
      {
        size_t (*myFileWrite)(int fileID, const void *restrict buffer, size_t len)
            = (size_t(*)(int, const void *restrict, size_t)) namespaceSwitchGet(NSSWITCH_FILE_WRITE).func;

        size_t nwrite = myFileWrite(fileID, gribbuffer, nbytes);
        if (nwrite != nbytes) SysError("Failed to write GRIB slice!");
      }
      break;
    }

  if (gribbuffer) Free(gribbuffer);
}

void
grb_write_var(stream_t *streamptr, int varID, int memtype, const void *data, size_t nmiss)
{
  int vlistID = streamptr->vlistID, gridID = vlistInqVarGrid(vlistID, varID), zaxisID = vlistInqVarZaxis(vlistID, varID),
      nlevs = zaxisInqSize(zaxisID);
  size_t gridsize = gridInqSize(gridID);
  double missval = vlistInqVarMissval(vlistID, varID);

  size_t chunkLen = gridsize;
  if (memtype == MEMTYPE_FLOAT)
    for (int levelID = 0; levelID < nlevs; levelID++)
      {
        const float *restrict fdata = ((const float *) data) + levelID * gridsize;

        size_t nmiss_slice = 0;
        if (nmiss)
          for (size_t i = 0; i < chunkLen; ++i) nmiss_slice += DBL_IS_EQUAL(fdata[i], missval);

        grb_write_var_slice(streamptr, varID, levelID, memtype, fdata, nmiss_slice);
      }
  else
    for (int levelID = 0; levelID < nlevs; levelID++)
      {
        const double *restrict ddata = ((const double *) data) + levelID * gridsize;

        size_t nmiss_slice = 0;
        if (nmiss)
          for (size_t i = 0; i < chunkLen; ++i) nmiss_slice += DBL_IS_EQUAL(ddata[i], missval);

        grb_write_var_slice(streamptr, varID, levelID, memtype, ddata, nmiss_slice);
      }
}

void
grb_write_record(stream_t *streamptr, int memtype, const void *data, size_t nmiss)
{
  int varID = streamptr->record->varID;
  int levelID = streamptr->record->levelID;
  grb_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss);
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
