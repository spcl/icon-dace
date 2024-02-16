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

// the single image implementation
static int
cdiStreamReadVar(int streamID, int varID, int memtype, void *data, size_t *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);
  const int filetype = streamptr->filetype;

  *nmiss = 0;

  if (memtype == MEMTYPE_FLOAT && cdiFiletypeIsExse(filetype)) return 1;

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_var(streamptr, varID, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvReadVarDP(streamptr, varID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extReadVarDP(streamptr, varID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegReadVarDP(streamptr, varID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_var(streamptr, varID, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype));
    }

  return status;
}

/*
@Function  streamReadVar
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, double *data, SizeType *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVar(int streamID, int varID, double *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  cdiStreamReadVar(streamID, varID, MEMTYPE_DOUBLE, data, &numMiss);
  *nmiss = (SizeType) numMiss;
}

/*
@Function  streamReadVarF
@Title     Read a variable

@Prototype void streamReadVar(int streamID, int varID, float *data, SizeType *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVar reads all the values of one time step of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarF(int streamID, int varID, float *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  if (cdiStreamReadVar(streamID, varID, MEMTYPE_FLOAT, data, &numMiss))
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      size_t elementCount = gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      elementCount *= (size_t) zaxisInqSize(vlistInqVarZaxis(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
      streamReadVar(streamID, varID, conversionBuffer, nmiss);
      for (size_t i = elementCount; i--;) data[i] = (float) conversionBuffer[i];
      Free(conversionBuffer);
    }
  *nmiss = (SizeType) numMiss;
}

static int
cdiStreamReadVarSlice(int streamID, int varID, int levelID, int memtype, void *data, size_t *nmiss)
{
  // May fail if memtype == MEMTYPE_FLOAT and the file format does not support single precision reading.
  // A value > 0 is returned in this case, otherwise it returns zero.
  int status = 0;

  if (CDI_Debug) Message("streamID = %d  varID = %d", streamID, varID);

  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);
  const int filetype = streamptr->filetype;

  *nmiss = 0;

  if (memtype == MEMTYPE_FLOAT && cdiFiletypeIsExse(filetype)) return 1;

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_var_slice(streamptr, varID, levelID, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srvReadVarSliceDP(streamptr, varID, levelID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: extReadVarSliceDP(streamptr, varID, levelID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: iegReadVarSliceDP(streamptr, varID, levelID, (double *) data, nmiss); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_var_slice(streamptr, varID, levelID, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(filetype));
    }

  return status;
}

/*
@Function  streamReadVarSlice
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSlice(int streamID, int varID, int levelID, double *data, SizeType *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSlice reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarSlice(int streamID, int varID, int levelID, double *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  if (cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_DOUBLE, data, &numMiss))
    {
      Warning("Unexpected error returned from cdiStreamReadVarSlice()!");
      size_t elementCount = gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      memset(data, 0, elementCount * sizeof(*data));
    }
  *nmiss = (SizeType) numMiss;
}

/*
@Function  streamReadVarSliceF
@Title     Read a horizontal slice of a variable

@Prototype void streamReadVarSliceF(int streamID, int varID, int levelID, float *data, SizeType *nmiss)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead}.
    @Item  varID     Variable identifier.
    @Item  levelID   Level identifier.
    @Item  data      Pointer to the location into which the data values are read.
                     The caller must allocate space for the returned values.
    @Item  nmiss     Number of missing values.

@Description
The function streamReadVarSliceF reads all the values of a horizontal slice of a variable
from an open dataset.
@EndFunction
*/
void
streamReadVarSliceF(int streamID, int varID, int levelID, float *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  if (cdiStreamReadVarSlice(streamID, varID, levelID, MEMTYPE_FLOAT, data, &numMiss))
    {
      // In case the file format does not support single precision reading,
      // we fall back to double precision reading, converting the data on the fly.
      size_t elementCount = gridInqSize(vlistInqVarGrid(streamInqVlist(streamID), varID));
      double *conversionBuffer = (double *) Malloc(elementCount * sizeof(*conversionBuffer));
      streamReadVarSlice(streamID, varID, levelID, conversionBuffer, nmiss);
      for (size_t i = elementCount; i--;) data[i] = (float) conversionBuffer[i];
      Free(conversionBuffer);
    }
  *nmiss = (SizeType) numMiss;
}

static void
stream_read_record(int streamID, int memtype, void *data, size_t *nmiss)
{
  check_parg(data);
  check_parg(nmiss);

  stream_t *streamptr = stream_to_pointer(streamID);

  *nmiss = 0;

  switch (cdiBaseFiletype(streamptr->filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB: grb_read_record(streamptr, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: srv_read_record(streamptr, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: ext_read_record(streamptr, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: ieg_read_record(streamptr, memtype, data, nmiss); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: cdf_read_record(streamptr, memtype, data, nmiss); break;
#endif
    default: Error("%s support not compiled in!", strfiletype(streamptr->filetype));
    }
}

void
streamReadRecord(int streamID, double *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  stream_read_record(streamID, MEMTYPE_DOUBLE, (void *) data, &numMiss);
  *nmiss = (SizeType) numMiss;
}

void
streamReadRecordF(int streamID, float *data, SizeType *nmiss)
{
  size_t numMiss = 0;
  stream_read_record(streamID, MEMTYPE_FLOAT, (void *) data, &numMiss);
  *nmiss = (SizeType) numMiss;
}
