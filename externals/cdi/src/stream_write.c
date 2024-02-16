#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "dmemory.h"
#include "namespace.h"

// the single image implementation
int
cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data, SizeType nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision writing.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0) Error("Writing of non-trivial subtypes not yet implemented!");

  // check taxis
  if (streamptr->curTsID == CDI_UNDEFID) streamDefTimestep(streamID, 0);

  const int filetype = streamptr->filetype;

  if (memtype == MEMTYPE_FLOAT && cdiFiletypeIsExse(filetype)) return 1;

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_write_var(streamptr, varID, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvWriteVarDP(streamptr, varID, (double *) data); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extWriteVarDP(streamptr, varID, (double *) data); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegWriteVarDP(streamptr, varID, (double *) data); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_write_var(streamptr, varID, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype));
    }

  return status;
}

/*
@Function  streamWriteVar
@Title     Write a variable

@Prototype void streamWriteVar(int streamID, int varID, const double *data, SizeType nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVar writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void
streamWriteVar(int streamID, int varID, const double *data, SizeType nmiss)
{
  void (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype, const void *data, SizeType nmiss)
      = (void (*)(int, int, int, const void *, SizeType)) namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;

  myCdiStreamWriteVar_(streamID, varID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

/*
@Function  streamWriteVarF
@Title     Write a variable

@Prototype void streamWriteVarF(int streamID, int varID, const float *data, SizeType nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarF writes the values of one time step of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void
streamWriteVarF(int streamID, int varID, const float *data, SizeType nmiss)
{
  int (*myCdiStreamWriteVar_)(int streamID, int varID, int memtype, const void *data, SizeType nmiss)
      = (int (*)(int, int, int, const void *, SizeType)) namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_).func;

  if (myCdiStreamWriteVar_(streamID, varID, MEMTYPE_FLOAT, (const void *) data, nmiss))
    {
      // In case the file format does not support single precision writing,
      // we fall back to double precision writing, converting the data on the fly.
      const int vlistID = streamInqVlist(streamID);
      SizeType elementCount = gridInqSize(vlistInqVarGrid(vlistID, varID));
      elementCount *= (SizeType) zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
      for (SizeType i = elementCount; i--;) conversionBuffer[i] = (double) data[i];
      myCdiStreamWriteVar_(streamID, varID, MEMTYPE_DOUBLE, (const void *) conversionBuffer, nmiss);
      Free(conversionBuffer);
    }
}

static int
cdiStreamWriteVarSlice(int streamID, int varID, int levelID, int memtype, const void *data, SizeType nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision writing.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d varID = %d", streamID, varID);

  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);
  if (subtypeInqActiveIndex(streamptr->vars[varID].subtypeID) != 0) Error("Writing of non-trivial subtypes not yet implemented!");

  // check taxis
  if (streamptr->curTsID == CDI_UNDEFID) streamDefTimestep(streamID, 0);

  const int filetype = streamptr->filetype;

  if (memtype == MEMTYPE_FLOAT && cdiFiletypeIsExse(filetype)) return 1;

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvWriteVarSliceDP(streamptr, varID, levelID, (double *) data); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extWriteVarSliceDP(streamptr, varID, levelID, (double *) data); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegWriteVarSliceDP(streamptr, varID, levelID, (double *) data); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_write_var_slice(streamptr, varID, levelID, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype));
    }

  return status;
}

/*
@Function  streamWriteVarSlice
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, SizeType nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSlice writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void
streamWriteVarSlice(int streamID, int varID, int levelID, const double *data, SizeType nmiss)
{
  cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

/*
@Function  streamWriteVarSliceF
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, SizeType nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to a block of single precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteVarSliceF writes the values of a horizontal slice of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void
streamWriteVarSliceF(int streamID, int varID, int levelID, const float *data, SizeType nmiss)
{
  if (cdiStreamWriteVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, (const void *) data, nmiss))
    {
      // In case the file format does not support single precision writing,
      // we fall back to double precision writing, converting the data on the fly.
      const SizeType elementCount = gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
      for (SizeType i = elementCount; i--;) conversionBuffer[i] = (double) data[i];
      streamWriteVarSlice(streamID, varID, levelID, conversionBuffer, nmiss);
      Free(conversionBuffer);
    }
}

void
streamWriteVarChunk(int streamID, int varID, const int rect[][2], const double *data, SizeType nmiss)
{
  void (*myCdiStreamWriteVarChunk_)(int streamID, int varID, int memtype, const int rect[3][2], const void *data, SizeType nmiss)
      = (void (*)(int, int, int, const int[3][2], const void *, SizeType)) namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_)
            .func;
  myCdiStreamWriteVarChunk_(streamID, varID, MEMTYPE_DOUBLE, rect, data, nmiss);
}

void
streamWriteVarChunkF(int streamID, int varID, const int rect[][2], const float *data, SizeType nmiss)
{
  void (*myCdiStreamWriteVarChunk_)(int streamID, int varID, int memtype, const int rect[3][2], const void *data, SizeType nmiss)
      = (void (*)(int, int, int, const int[3][2], const void *, SizeType)) namespaceSwitchGet(NSSWITCH_STREAM_WRITE_VAR_CHUNK_)
            .func;
  myCdiStreamWriteVarChunk_(streamID, varID, MEMTYPE_FLOAT, rect, data, nmiss);
}

// single image implementation
void
cdiStreamWriteVarChunk_(int streamID, int varID, int memtype, const int rect[][2], const void *data, SizeType nmiss)
{
  if (CDI_Debug) Message("streamID = %d varID = %d", streamID, varID);

  stream_t *streamptr = stream_to_pointer(streamID);

  // streamDefineTaxis(streamID);

  const int filetype = streamptr->filetype;

  switch (cdiBaseFiletype(filetype))
    {
#if defined(HAVE_LIBGRIB)
    case CDI_FILETYPE_GRIB:
#endif
#if defined(HAVE_LIBSERVICE)
    case CDI_FILETYPE_SRV:
#endif
#if defined(HAVE_LIBEXTRA)
    case CDI_FILETYPE_EXT:
#endif
#if defined(HAVE_LIBIEG)
    case CDI_FILETYPE_IEG:
#endif
#if defined(HAVE_LIBGRIB) || defined(HAVE_LIBSERVICE) || defined(HAVE_LIBEXTRA) || defined(HAVE_LIBIEG)
      xabort("streamWriteVarChunk not implemented for filetype %s!", strfiletype(filetype));
      break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_write_var_chunk(streamptr, varID, memtype, rect, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype)); break;
    }
}

static void
stream_write_record(int streamID, int memtype, const void *data, SizeType nmiss)
{
  check_parg(data);

  stream_t *streamptr = stream_to_pointer(streamID);

  switch (cdiBaseFiletype(streamptr->filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_write_record(streamptr, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srv_write_record(streamptr, memtype, data); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: ext_write_record(streamptr, memtype, data); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: ieg_write_record(streamptr, memtype, data); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_write_record(streamptr, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(streamptr->filetype));
    }
}

/*
@Function  streamWriteRecord
@Title     Write a horizontal slice of a variable

@Prototype void streamWriteRecord(int streamID, const double *data, SizeType nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  data      Pointer to a block of double precision floating point data values to be written.
    @Item  nmiss     Number of missing values.

@Description
The function streamWriteRecord writes the values of a horizontal slice (record) of a variable to an open dataset.
The values are converted to the external data type of the variable, if necessary.
@EndFunction
*/
void
streamWriteRecord(int streamID, const double *data, SizeType nmiss)
{
  stream_write_record(streamID, MEMTYPE_DOUBLE, (const void *) data, nmiss);
}

void
streamWriteRecordF(int streamID, const float *data, SizeType nmiss)
{
  stream_write_record(streamID, MEMTYPE_FLOAT, (const void *) data, nmiss);
}
