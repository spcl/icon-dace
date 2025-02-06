#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#ifdef HAVE_LIBFDB5
#include "cdi_fdb.h"
#endif

#include <sys/stat.h>  // struct stat
#include <ctype.h>
#include <string.h>

#include "binary.h"
#include "cdi.h"
#include "cdi_across.h"
#include "cdi_int.h"
#include "cdi_cksum.h"
#include "cdf.h"
#include "dmemory.h"
#include "error.h"
#include "stream_cgribex.h"
#include "stream_grb.h"
#include "stream_cdf.h"
#include "stream_srv.h"
#include "stream_ext.h"
#include "stream_ieg.h"
#include "file.h"
#include "cgribex.h"
#include "gribapi.h"
#include "vlist.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "async_worker.h"

static stream_t *stream_new_entry(int resH);
static int streamCompareP(void *streamptr1, void *streamptr2);
static void streamDestroyP(void *streamptr);
static void streamPrintP(void *streamptr, FILE *fp);
static int streamGetPackSize(void *streamptr, void *context);
static void streamPack(void *streamptr, void *buff, int size, int *position, void *context);
static int streamTxCode(void *streamptr);

const resOps streamOps = { streamCompareP, streamDestroyP, streamPrintP, streamGetPackSize, streamPack, streamTxCode };

static int
getByteorder(int byteswap)
{
  int byteorder = -1;

  switch (HOST_ENDIANNESS)
    {
    case CDI_BIGENDIAN: byteorder = byteswap ? CDI_LITTLEENDIAN : CDI_BIGENDIAN; break;
    case CDI_LITTLEENDIAN: byteorder = byteswap ? CDI_BIGENDIAN : CDI_LITTLEENDIAN; break;
    /* FIXME: does not currently adjust for PDP endianness */
    case CDI_PDPENDIAN:
    default: Error("unhandled endianness");
    }

  return byteorder;
}

// used also in CDO
int
cdiGetFiletype(const char *uri, int *byteorder)
{
  // clang-format off
  int filetype = CDI_EUFTYPE;
  int swap = 0;
  int version;
  long recpos;

  const char *filename;
  int protocol = cdiGetProtocol(uri, &filename);

  switch (protocol)
    {
    case CDI_PROTOCOL_ACROSS: return CDI_FILETYPE_GRB2;
    case CDI_PROTOCOL_FDB:    return CDI_FILETYPE_GRB2;
    case CDI_PROTOCOL_OTHER:  return CDI_FILETYPE_NC4;  // support for NetCDF remote types and ESDM
    case CDI_PROTOCOL_FILE:
      // handled below;
      break;
    }

  int fileID = fileOpen(filename, "r");

  if      (fileID == CDI_UNDEFID) return CDI_ESYSTEM;
  else if (fileID == -2)          return CDI_ETMOF;

  char buffer[8];
  if (fileRead(fileID, buffer, 8) != 8)
    {
      struct stat buf;
      if (stat(filename, &buf) == 0)
        {
          if (buf.st_size == 0)      return CDI_EISEMPTY;
          if (buf.st_mode & S_IFDIR) return CDI_EISDIR;
        }

      return CDI_EUFTYPE;
    }

  fileRewind(fileID);

  if (memcmp(buffer, "GRIB", 4) == 0)
    {
      version = buffer[7];
      if      (version <= 1) filetype = CDI_FILETYPE_GRB;
      else if (version == 2) filetype = CDI_FILETYPE_GRB2;
    }
  else if (memcmp(buffer, "CDF\001", 4) == 0) { filetype = CDI_FILETYPE_NC; }
  else if (memcmp(buffer, "CDF\002", 4) == 0) { filetype = CDI_FILETYPE_NC2; }
  else if (memcmp(buffer, "CDF\005", 4) == 0) { filetype = CDI_FILETYPE_NC5; }
  else if (memcmp(buffer + 1, "HDF", 3) == 0) { filetype = CDI_FILETYPE_NC4; }
#ifdef HAVE_LIBSERVICE
  else if (srvCheckFiletype(fileID, &swap))   { filetype = CDI_FILETYPE_SRV; }
#endif
#ifdef HAVE_LIBEXTRA
  else if (extCheckFiletype(fileID, &swap))   { filetype = CDI_FILETYPE_EXT; }
#endif
#ifdef HAVE_LIBIEG
  else if (iegCheckFiletype(fileID, &swap))   { filetype = CDI_FILETYPE_IEG; }
#endif
#ifdef HAVE_LIBGRIB
  else if (gribCheckSeek(fileID, &recpos, &version) == 0)
    {
      if      (version <= 1) filetype = CDI_FILETYPE_GRB;
      else if (version == 2) filetype = CDI_FILETYPE_GRB2;
    }
#endif
  // clang-format on

  if (CDI_Debug && filetype != CDI_EUFTYPE) Message("found %s file = %s", strfiletype(filetype), filename);

  fileClose(fileID);

  *byteorder = getByteorder(swap);

  return filetype;
}

/*
@Function  streamInqFiletype
@Title     Get the filetype

@Prototype int streamInqFiletype(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqFiletype} returns the filetype of a stream.

@Result
@func{streamInqFiletype} returns the type of the file format,
one of the set of predefined CDI file format types.
The valid CDI file format types are @func{CDI_FILETYPE_GRB}, @func{CDI_FILETYPE_GRB2}, @func{CDI_FILETYPE_NC},
@func{CDI_FILETYPE_NC2}, @func{CDI_FILETYPE_NC4}, @func{CDI_FILETYPE_NC4C}, @func{CDI_FILETYPE_NC5},
@func{CDI_FILETYPE_NCZARR}, @func{CDI_FILETYPE_SRV}, @func{CDI_FILETYPE_EXT} and @func{CDI_FILETYPE_IEG}.
@EndFunction
*/
int
streamInqFiletype(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->filetype;
}

int
getByteswap(int byteorder)
{
  int byteswap = -1;

  switch (byteorder)
    {
    case CDI_BIGENDIAN:
    case CDI_LITTLEENDIAN:
    case CDI_PDPENDIAN: byteswap = (HOST_ENDIANNESS != byteorder); break;
    case -1: break;
    default: Error("unexpected byteorder %d query!", byteorder);
    }

  return byteswap;
}

void
streamDefMaxSteps(int streamID, int maxSteps)
{
  if (maxSteps >= 0)
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->maxSteps = maxSteps;
    }
}

int
streamInqNumSteps(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  long ntsteps = streamptr->ntsteps;
  if (ntsteps == (long) CDI_UNDEFID)
    {
      int tsID = 0;
      while (streamInqTimestep(streamID, tsID++)) ntsteps = streamptr->ntsteps;
    }

  return (int) ntsteps;
}

static long
get_max_global_recs(stream_t *streamptr)
{
  long maxGlobalRecs = -1;
  const tsteps_t *tsteps = streamptr->tsteps;
  if (tsteps)
    {
      maxGlobalRecs = tsteps[0].nrecs;
      long numSteps = streamptr->ntsteps;
      if (numSteps > 1) maxGlobalRecs += tsteps[1].nrecs * (numSteps - 1);
    }
  return maxGlobalRecs;
}

void
streamDefNumWorker(int streamID, int numWorker)
{
  if (numWorker > 0)
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      int filetype = streamptr->filetype;

      if (streamptr->filemode == 'r')
        {
          if (cdiBaseFiletype(filetype) == CDI_FILETYPE_GRIB)
            {
              (void) streamInqNumSteps(streamID);
              streamptr->maxGlobalRecs = get_max_global_recs(streamptr);
            }
#ifdef HAVE_LIBNETCDF
          else if (filetype == CDI_FILETYPE_NCZARR || (CDI_Test && cdiBaseFiletype(filetype) == CDI_FILETYPE_NETCDF))
            {
              streamptr->maxGlobalRecs = get_max_global_recs(streamptr);
              if (CDI_Test) Message("numWorker=%d", numWorker);
              if (CDI_Test) Message("maxGlobalRecs=%ld", streamptr->maxGlobalRecs);
              if (streamptr->maxGlobalRecs == -1) xabort("Internal error: number of timesteps missing!");
              if (streamptr->maxGlobalRecs == 1) numWorker = 0;
              if (numWorker > streamptr->maxGlobalRecs) numWorker = streamptr->maxGlobalRecs;
              if (streamptr->chunkSizeTdim > 1 && numWorker > streamptr->nvars) numWorker = streamptr->nvars;
              if (streamptr->chunkSizeZdim > 1) numWorker = 0;
              if (CDI_Test) Message("chunkSizeTdim=%d chunkSizeZdim=%d", streamptr->chunkSizeTdim, streamptr->chunkSizeZdim);
            }
#endif
          else
            {
              numWorker = 0;
            }

          streamptr->numWorker = numWorker;
          if (CDI_Debug || CDI_Test) Message("Number of asynchronous worker: %d", numWorker);
        }
    }
}

/*
@Function  streamDefByteorder
@Title     Define the byte order

@Prototype void streamDefByteorder(int streamID, int byteorder)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  byteorder The byte order of a dataset, one of the CDI constants @func{CDI_BIGENDIAN} and
                     @func{CDI_LITTLEENDIAN}.

@Description
The function @func{streamDefByteorder} defines the byte order of a binary dataset
with the file format type @func{CDI_FILETYPE_SRV}, @func{CDI_FILETYPE_EXT} or @func{CDI_FILETYPE_IEG}.

@EndFunction
*/
void
streamDefByteorder(int streamID, int byteorder)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  streamptr->byteorder = byteorder;
  int filetype = streamptr->filetype;

  switch (filetype)
    {
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV:
      {
        srvrec_t *srvp = (srvrec_t *) streamptr->record->objectp;
        srvp->byteswap = getByteswap(byteorder);

        break;
      }
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT:
      {
        extrec_t *extp = (extrec_t *) streamptr->record->objectp;
        extp->byteswap = getByteswap(byteorder);

        break;
      }
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG:
      {
        iegrec_t *iegp = (iegrec_t *) streamptr->record->objectp;
        iegp->byteswap = getByteswap(byteorder);

        break;
      }
#endif
    }

  reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
}

/*
@Function  streamInqByteorder
@Title     Get the byte order

@Prototype int streamInqByteorder(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqByteorder} returns the byte order of a binary dataset
with the file format type @func{CDI_FILETYPE_SRV}, @func{CDI_FILETYPE_EXT} or @func{CDI_FILETYPE_IEG}.

@Result
@func{streamInqByteorder} returns the type of the byte order.
The valid CDI byte order types are @func{CDI_BIGENDIAN} and @func{CDI_LITTLEENDIAN}

@EndFunction
*/
int
streamInqByteorder(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->byteorder;
}

const char *
streamFilesuffix(int filetype)
{
  static const char *noSuffix = "";
  static const char *ncSuffix = ".nc";
  static const char *grbSuffix = ".grb";
  static const char *srvSuffix = ".srv";
  static const char *extSuffix = ".ext";
  static const char *iegSuffix = ".ieg";

  switch (cdiBaseFiletype(filetype))
    {
    case CDI_FILETYPE_GRIB: return grbSuffix;
    case CDI_FILETYPE_SRV: return srvSuffix;
    case CDI_FILETYPE_EXT: return extSuffix;
    case CDI_FILETYPE_IEG: return iegSuffix;
    case CDI_FILETYPE_NETCDF: return ncSuffix;
    default: return noSuffix;
    }
}

const char *
streamFilename(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->filename;
}

static int
cdiInqContents(stream_t *streamptr)
{
  int status = 0;
  int filetype = streamptr->filetype;

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB:
      {
        switch (streamptr->protocol)
          {
          case CDI_PROTOCOL_FDB: status = fdbInqContents(streamptr); break;

          case CDI_PROTOCOL_ACROSS:  // TODO read from ACROSS
          case CDI_PROTOCOL_OTHER:
          case CDI_PROTOCOL_FILE: status = grbInqContents(streamptr); break;
          }
        break;
      }
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV: status = srvInqContents(streamptr); break;
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT: status = extInqContents(streamptr); break;
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG: status = iegInqContents(streamptr); break;
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF: status = cdfInqContents(streamptr); break;
#endif
    default:
      {
        status = CDI_ELIBNAVAIL;
        if (CDI_Debug) Message("%s support not compiled in!", strfiletype(filetype));
      }
    }

  if (status == 0)
    {
      int taxisID = vlistInqTaxis(streamptr->vlistID);
      if (taxisID != CDI_UNDEFID)
        {
          taxis_t *taxisptr1 = &streamptr->tsteps[0].taxis;
          taxis_t *taxisptr2 = taxisPtr(taxisID);
          ptaxisCopy(taxisptr2, taxisptr1);
        }
    }

  return status;
}

int
cdiGetProtocol(const char *uri, const char **filename)
{
  const char *pos = strstr(uri, "://");
  if (pos == NULL)
    {
      *filename = uri;
      return CDI_PROTOCOL_FILE;
    }

  int protocollen = pos - uri;
  *filename = pos + 3;

  // if (strncmp(uri, "file", protocollen) == 0) return CDI_PROTOCOL_FILE; // file is already used in NetCDF
  if (strncmp(uri, "fdb", protocollen) == 0) return CDI_PROTOCOL_FDB;
  if (strncmp(uri, "across", protocollen) == 0) return CDI_PROTOCOL_ACROSS;

  *filename = uri;

  return CDI_PROTOCOL_OTHER;
}

int
cdiStreamOpenDefaultDelegate(const char *uri, char filemode, int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
{
  int fileID;
  char temp[2] = { filemode, 0 };

  const char *filename;
  streamptr->protocol = cdiGetProtocol(uri, &filename);

  switch (streamptr->protocol)
    {
    case CDI_PROTOCOL_ACROSS:
#if defined(HAVE_ACROSS) && defined(HAVE_LIBGRIB_API)
      if (filetype != CDI_FILETYPE_GRB2)
        {
          Warning("ACROSS needs to be used with GRIB2");
          return CDI_EUFTYPE;
        }
      fileID = across_connect(filename, filemode, streamptr);
      if (fileID >= 0)
        {
          streamptr->filetype = filetype;
          if (recordBufIsToBeCreated)
            {
              streamptr->record = (Record *) Malloc(sizeof(Record));
              streamptr->record->buffer = NULL;
            }
        }
      return fileID;
#else
#ifdef HAVE_ACROSS
      Warning("ecCodes support not compiled in (Needed for ACROSS)!");
#else
      Warning("ACROSS support not compiled in!");
#endif
      return CDI_ELIBNAVAIL;
#endif

    case CDI_PROTOCOL_FDB:
#if defined(HAVE_LIBFDB5) && defined(HAVE_LIBGRIB_API)

      if (filetype != CDI_FILETYPE_GRB && filetype != CDI_FILETYPE_GRB2)
        {
          Warning("FDB5 needs to be used with GRIB or GRIB2");
          return CDI_EUFTYPE;
        }

      check_fdb_error(fdb_initialise());
      check_fdb_error(fdb_new_handle((fdb_handle_t **) &(streamptr->protocolData)));
      streamptr->filetype = filetype;
      if (recordBufIsToBeCreated)
        {
          streamptr->record = (Record *) Malloc(sizeof(Record));
          streamptr->record->buffer = NULL;
        }
      return 88;

#else  // !(defined(HAVE_LIBFDB5) && defined(HAVE_LIBGRIB_API))

#ifdef HAVE_LIBFDB5
      Warning("ecCodes support not compiled in (Needed for FDB5)!");
#else
      Warning("FDB5 support not compiled in!");
#endif

      return CDI_ELIBNAVAIL;
#endif

    case CDI_PROTOCOL_OTHER:
    case CDI_PROTOCOL_FILE:
      // handled below;
      break;
    }

  switch (filetype)
    {
#if defined(HAVE_LIBGRIB) && (defined(HAVE_LIBCGRIBEX) || defined(HAVE_LIBGRIB_API))
    case CDI_FILETYPE_GRB:
      {
        fileID = gribOpen(filename, temp);
        if (fileID < 0) return CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) Malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
#ifdef HAVE_LIBCGRIBEX
            streamptr->record->objectp = cgribexNew();
#else
            streamptr->record->objectp = NULL;
#endif
          }
        break;
      }
#ifdef HAVE_LIBGRIB_API
    case CDI_FILETYPE_GRB2:
      {
        fileID = gribOpen(filename, temp);
        if (fileID < 0) return CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) Malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
          }
        break;
      }
#endif
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV:
      {
        fileID = fileOpen(filename, temp);
        if (fileID < 0) return CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) Malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->objectp = srvNew();
          }
        break;
      }
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT:
      {
        fileID = fileOpen(filename, temp);
        if (fileID < 0) return CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) Malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->objectp = extNew();
          }
        break;
      }
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG:
      {
        fileID = fileOpen(filename, temp);
        if (fileID < 0) return CDI_ESYSTEM;
        if (recordBufIsToBeCreated)
          {
            streamptr->record = (Record *) Malloc(sizeof(Record));
            streamptr->record->buffer = NULL;
            streamptr->record->objectp = iegNew();
          }
        break;
      }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC5:
      {
        fileID = cdfOpen(filename, temp, filetype);
        break;
      }
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
    case CDI_FILETYPE_NCZARR:
      {
        fileID = cdf4Open(filename, temp, &filetype);
        break;
      }
#endif
    default:
      {
        if (CDI_Debug) Message("%s support not compiled in!", strfiletype(filetype));
        return CDI_ELIBNAVAIL;
      }
    }

  streamptr->filetype = filetype;

  return fileID;
}

static int
stream_create_vlist(stream_t *streamptr, CdiQuery *query)
{
  int vlistID = vlistCreate();
  if (vlistID < 0) return CDI_ELIMIT;

  cdiVlistMakeInternal(vlistID);
  streamptr->vlistID = vlistID;

  if (query) streamptr->query = query;

  int status = cdiInqContents(streamptr);
  if (status >= 0)
    {
      vlist_t *vlistptr = vlist_to_pointer(streamptr->vlistID);
      vlistptr->ntsteps = streamptr->ntsteps;
      cdiVlistMakeImmutable(vlistID);
    }

  return status;
}

int
streamOpenID(const char *filename, char filemode, int filetype, int resH)
{
  if (CDI_Debug) Message("Open %s mode %c file %s", strfiletype(filetype), filemode, filename ? filename : "(NUL)");

  if (!filename || filetype < 0) return CDI_EINVAL;

  stream_t *streamptr = stream_new_entry(resH);
  int streamID = CDI_ESYSTEM;

#ifndef HAVE_NC4HDF5_THREADSAFE
  if (CDI_Threadsafe)
    {
#ifndef HAVE_LIBPTHREAD
      Error("CDI threadsafe failed, pthread support not compiled!");
#endif
      if (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C) streamptr->lockIO = true;
    }
#endif

  if (streamptr->lockIO) CDI_IO_LOCK();

  int (*streamOpenDelegate)(const char *filename, char filemode, int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
      = (int (*)(const char *, char, int, stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_OPEN_BACKEND).func;

  int fileID = streamOpenDelegate(filename, filemode, filetype, streamptr, 1);
  if (fileID < 0)
    {
      streamID = fileID;
      if (streamptr->record) Free(streamptr->record);
      reshRemove(streamptr->self, &streamOps);
      Free(streamptr);
    }
  else
    {
      streamID = streamptr->self;
      if (streamID < 0) return CDI_ELIMIT;

      streamptr->filemode = filemode;
      streamptr->filename = strdup(filename);
      streamptr->fileID = fileID;
    }

  if (streamptr->lockIO) CDI_IO_UNLOCK();

  return streamID;
}

static int
streamOpen(const char *filename, const char *filemode, int filetype)
{
  if (!filemode || strlen(filemode) != 1) return CDI_EINVAL;
  return streamOpenID(filename, (char) tolower(filemode[0]), filetype, CDI_UNDEFID);
}

static int
streamOpenA(const char *filename, const char *filemode, int filetype)
{
  if (CDI_Debug) Message("Open %s file (mode=%c); filename: %s", strfiletype(filetype), (int) *filemode, filename);
  if (CDI_Debug) printf("streamOpenA: %s\n", filename);  // seg fault without this line on thunder/squall with "cdo cat x y"

  if (!filename || !filemode || filetype < 0) return CDI_EINVAL;

  stream_t *streamptr = stream_new_entry(CDI_UNDEFID);
  int fileID = CDI_UNDEFID;

  {
    int (*streamOpenDelegate)(const char *filename, char filemode, int filetype, stream_t *streamptr, int recordBufIsToBeCreated)
        = (int (*)(const char *, char, int, stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_OPEN_BACKEND).func;

    fileID = streamOpenDelegate(filename, 'r', filetype, streamptr, 1);
  }

  if (fileID == CDI_UNDEFID || fileID == CDI_ELIBNAVAIL || fileID == CDI_ESYSTEM) return fileID;

  int streamID = streamptr->self;

  streamptr->filemode = tolower(*filemode);
  streamptr->filename = strdup(filename);
  streamptr->fileID = fileID;

  streamptr->vlistID = vlistCreate();
  cdiVlistMakeInternal(streamptr->vlistID);
  // cdiReadByteorder(streamID);
  int status = cdiInqContents(streamptr);
  if (status < 0) return status;
  vlist_t *vlistptr = vlist_to_pointer(streamptr->vlistID);
  vlistptr->ntsteps = streamInqNumSteps(streamID);

  // Needed for NetCDF4
  for (int varID = 0; varID < vlistptr->nvars; ++varID) streamptr->vars[varID].defmiss = true;

  if (str_is_equal(filemode, "r")) cdiVlistMakeImmutable(streamptr->vlistID);

  {
    void (*streamCloseDelegate)(stream_t * streamptr, int recordBufIsToBeDeleted)
        = (void (*)(stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_CLOSE_BACKEND).func;

    streamCloseDelegate(streamptr, 0);
  }

  switch (filetype)
    {
#if defined(HAVE_LIBGRIB) && (defined(HAVE_LIBCGRIBEX) || defined(HAVE_LIBGRIB_API))
    case CDI_FILETYPE_GRB:
#ifdef HAVE_LIBGRIB_API
    case CDI_FILETYPE_GRB2:
#endif
      {
        fileID = gribOpen(filename, filemode);
        if (fileID != CDI_UNDEFID) gribContainersNew(streamptr);
        break;
      }
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV:
      {
        fileID = fileOpen(filename, filemode);
        break;
      }
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT:
      {
        fileID = fileOpen(filename, filemode);
        break;
      }
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG:
      {
        fileID = fileOpen(filename, filemode);
        break;
      }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NC:
    case CDI_FILETYPE_NC2:
    case CDI_FILETYPE_NC5:
      {
        fileID = cdfOpen(filename, filemode, filetype);
        streamptr->ncmode = 2;
        break;
      }
    case CDI_FILETYPE_NC4:
    case CDI_FILETYPE_NC4C:
      {
        fileID = cdf4Open(filename, filemode, &filetype);
        streamptr->ncmode = 2;
        break;
      }
    case CDI_FILETYPE_NCZARR:
      {
        Warning("%s not available in append mode!", strfiletype(filetype));
        return CDI_ELIBNAVAIL;
      }
#endif
    default:
      {
        if (CDI_Debug) Message("%s support not compiled in!", strfiletype(filetype));
        return CDI_ELIBNAVAIL;
      }
    }

  if (fileID == CDI_UNDEFID)
    streamID = CDI_UNDEFID;
  else
    streamptr->fileID = fileID;

  return streamID;
}

/*
@Function  streamOpenRead
@Title     Open a dataset for reading

@Prototype int streamOpenRead(const char *path)
@Parameter
    @Item  path  The name of the dataset to be read.

@Description
The function @func{streamOpenRead} opens an existing dataset for reading.

@Result
Upon successful completion @func{streamOpenRead} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenRead} to open an existing NetCDF
file named @func{foo.nc} for reading:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenRead("foo.nc");
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int
streamOpenRead(const char *filename)
{
  cdiInitialize();

  int byteorder = 0;
  int filetype = cdiGetFiletype(filename, &byteorder);
  if (filetype < 0) return filetype;

  int streamID = streamOpen(filename, "r", filetype);
  if (streamID >= 0)
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;

      int status = stream_create_vlist(streamptr, NULL);
      if (status < 0)
        {
          streamID = status;
          if (streamptr->record) Free(streamptr->record);
          reshRemove(streamptr->self, &streamOps);
        }
    }

  return streamID;
}

int
streamOpenReadQuery(const char *filename, CdiQuery *query)
{
  cdiInitialize();

  int byteorder = 0;
  int filetype = cdiGetFiletype(filename, &byteorder);
  if (filetype < 0) return filetype;

  if (cdiBaseFiletype(filetype) != CDI_FILETYPE_NETCDF && filetype != CDI_FILETYPE_GRB2) return CDI_EQNAVAIL;

  int streamID = streamOpen(filename, "r", filetype);
  if (streamID >= 0)
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;

      int status = stream_create_vlist(streamptr, query);
      if (status < 0)
        {
          streamID = status;
          if (streamptr->record) Free(streamptr->record);
          reshRemove(streamptr->self, &streamOps);
        }
    }

  return streamID;
}

int
streamOpenAppend(const char *filename)
{
  cdiInitialize();

  int byteorder = 0;
  int filetype = cdiGetFiletype(filename, &byteorder);
  if (filetype < 0) return filetype;

  int streamID = streamOpenA(filename, "a", filetype);
  if (streamID >= 0)
    {
      stream_t *streamptr = stream_to_pointer(streamID);
      streamptr->byteorder = byteorder;
    }

  return streamID;
}

/*
@Function  streamOpenWrite
@Title     Create a new dataset

@Prototype int streamOpenWrite(const char *path, int filetype)
@Parameter
    @Item  path      The name of the new dataset.
    @Item  filetype  The type of the file format, one of the set of predefined CDI file format types.
                     The valid CDI file format types are @func{CDI_FILETYPE_GRB}, @func{CDI_FILETYPE_GRB2}, @func{CDI_FILETYPE_NC},
                     @func{CDI_FILETYPE_NC2}, @func{CDI_FILETYPE_NC4}, @func{CDI_FILETYPE_NC4C}, @func{CDI_FILETYPE_NC5},
                     @func{CDI_FILETYPE_NCZARR}, @func{CDI_FILETYPE_SRV}, @func{CDI_FILETYPE_EXT} and @func{CDI_FILETYPE_IEG}.

@Description
The function @func{streamOpenWrite} creates a new datset.
@Result
Upon successful completion @func{streamOpenWrite} returns an identifier to the
open stream. Otherwise, a negative number with the error status is returned.

@Errors
@List
   @Item  CDI_ESYSTEM     Operating system error.
   @Item  CDI_EINVAL      Invalid argument.
   @Item  CDI_EUFILETYPE  Unsupported file type.
   @Item  CDI_ELIBNAVAIL  Library support not compiled in.
@EndList

@Example
Here is an example using @func{streamOpenWrite} to create a new NetCDF file named @func{foo.nc} for writing:

@Source
#include "cdi.h"
   ...
int streamID;
   ...
streamID = streamOpenWrite("foo.nc", CDI_FILETYPE_NC);
if ( streamID < 0 ) handle_error(streamID);
   ...
@EndSource
@EndFunction
*/
int
streamOpenWrite(const char *filename, int filetype)
{
  cdiInitialize();

  return streamOpen(filename, "w", filetype);
}

static void
streamDefaultValue(stream_t *streamptr)
{
  streamptr->self = CDI_UNDEFID;
  streamptr->accesstype = CDI_UNDEFID;
  streamptr->accessmode = 0;
  streamptr->filetype = CDI_FILETYPE_UNDEF;
  streamptr->byteorder = CDI_UNDEFID;
  streamptr->fileID = 0;
  streamptr->filemode = 0;
  streamptr->numvals = 0;
  streamptr->filename = NULL;
  streamptr->record = NULL;
  streamptr->query = NULL;
  streamptr->varsAllocated = 0;
  streamptr->nrecs = 0;
  streamptr->nvars = 0;
  streamptr->vars = NULL;
  streamptr->ncmode = 0;
  streamptr->curTsID = CDI_UNDEFID;
  streamptr->rtsteps = 0;
  streamptr->ntsteps = CDI_UNDEFID;
  streamptr->maxSteps = CDI_UNDEFID;
  streamptr->tsteps = NULL;
  streamptr->tstepsTableSize = 0;
  streamptr->tstepsNextID = 0;
  streamptr->vlistID = CDI_UNDEFID;
  streamptr->globalatts = 0;
  streamptr->localatts = 0;
  streamptr->unreduced = cdiDataUnreduced;
  streamptr->have_missval = cdiHaveMissval;
  streamptr->comptype = CDI_COMPRESS_NONE;
  streamptr->complevel = 0;
  streamptr->shuffle = 0;
  streamptr->sortname = (cdiSortName > 0);
  streamptr->lockIO = CDI_Lock_IO;
  // netcdf4/HDF5 filter
  streamptr->filterId = 0;
  streamptr->numParams = 0;
  streamptr->maxParams = sizeof(streamptr->params) / sizeof(streamptr->params[0]);

  basetimeInit(&streamptr->basetime);

#ifdef HAVE_LIBNETCDF
  streamptr->nc_complex_float_id = CDI_UNDEFID;
  streamptr->nc_complex_double_id = CDI_UNDEFID;

  for (int i = 0; i < MAX_ZAXES_PS; i++) streamptr->zaxisID[i] = CDI_UNDEFID;
  for (int i = 0; i < MAX_ZAXES_PS; i++) streamptr->nczvarID[i] = CDI_UNDEFID;

  for (int i = 0; i < MAX_GRIDS_PS; i++)
    {
      streamptr->ncgrid[i].start = CDI_UNDEFID;
      streamptr->ncgrid[i].count = CDI_UNDEFID;
      streamptr->ncgrid[i].gridID = CDI_UNDEFID;
      for (size_t j = 0; j < CDF_SIZE_ncIDs; ++j) streamptr->ncgrid[i].ncIDs[j] = CDI_UNDEFID;
    }

  streamptr->ncNumDims = 0;
  for (int i = 0; i < MAX_DIMS_PS; i++) streamptr->ncDimID[i] = CDI_UNDEFID;
  for (int i = 0; i < MAX_DIMS_PS; i++) streamptr->ncDimLen[i] = 0;

  streamptr->vct.ilev = 0;
  streamptr->vct.mlev = 0;
  streamptr->vct.ilevID = CDI_UNDEFID;
  streamptr->vct.mlevID = CDI_UNDEFID;

  streamptr->chunkSizeTdim = 0;
  streamptr->chunkSizeZdim = 0;
#endif
  streamptr->maxGlobalRecs = CDI_UNDEFID;

  streamptr->gribContainers = NULL;

  streamptr->numWorker = 0;
  streamptr->nextGlobalRecId = 0;
  streamptr->cachedTsID = -1;
  streamptr->jobs = NULL;
  streamptr->jobManager = NULL;

  streamptr->protocolData = NULL;

#ifdef HAVE_LIBFDB5
  streamptr->fdbNumItems = 0;
  streamptr->fdbKeyValueList = NULL;
#endif
}

static stream_t *
stream_new_entry(int resH)
{
  cdiInitialize(); /* ***************** make MT version !!! */

  stream_t *streamptr = (stream_t *) Malloc(sizeof(stream_t));
  streamDefaultValue(streamptr);

  if (resH == CDI_UNDEFID)
    streamptr->self = reshPut(streamptr, &streamOps);
  else
    {
      streamptr->self = resH;
      reshReplace(resH, streamptr, &streamOps);
    }

  return streamptr;
}

void
cdiStreamCloseDefaultDelegate(stream_t *streamptr, int recordBufIsToBeDeleted)
{
  int fileID = streamptr->fileID;
  int filetype = streamptr->filetype;

  switch (streamptr->protocol)
    {
    case CDI_PROTOCOL_ACROSS:
#ifdef HAVE_ACROSS
      if (fileID) across_disconnect(fileID);
      if (streamptr->protocolData)
        {
          Free(((across_info_t *) streamptr->protocolData)->expid);
          Free(streamptr->protocolData);
          streamptr->protocolData = NULL;
        }
#endif
      return;

    case CDI_PROTOCOL_FDB:
#ifdef HAVE_LIBFDB5
      if (streamptr->protocolData) check_fdb_error(fdb_delete_handle(streamptr->protocolData));
      streamptr->protocolData = NULL;
#endif
      return;

    case CDI_PROTOCOL_OTHER:
    case CDI_PROTOCOL_FILE:
      // handled below;
      break;
    }

  if (fileID == CDI_UNDEFID)
    {
      Warning("File %s not open!", streamptr->filename);
      return;
    }

  switch (cdiBaseFiletype(filetype))
    {
#if defined(HAVE_LIBGRIB) && (defined(HAVE_LIBCGRIBEX) || defined(HAVE_LIBGRIB_API))
    case CDI_FILETYPE_GRIB:
      if (filetype == CDI_FILETYPE_GRB)
        {
          gribClose(fileID);
          if (recordBufIsToBeDeleted) gribContainersDelete(streamptr);
#ifdef HAVE_LIBCGRIBEX
          if (recordBufIsToBeDeleted) cgribexDelete(streamptr->record->objectp);
#endif
        }
      else if (filetype == CDI_FILETYPE_GRB2)
        {
          gribClose(fileID);
          if (recordBufIsToBeDeleted) gribContainersDelete(streamptr);
        }
      break;
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV:
      {
        fileClose(fileID);
        if (recordBufIsToBeDeleted) srvDelete(streamptr->record->objectp);
        break;
      }
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT:
      {
        fileClose(fileID);
        if (recordBufIsToBeDeleted) extDelete(streamptr->record->objectp);
        break;
      }
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG:
      {
        fileClose(fileID);
        if (recordBufIsToBeDeleted) iegDelete(streamptr->record->objectp);
        break;
      }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF:
      {
        cdfClose(fileID);
        if (streamptr->ntsteps == 0 && streamptr->tsteps != NULL)
          {
            if (streamptr->tsteps[0].records)
              {
                Free(streamptr->tsteps[0].records);
                streamptr->tsteps[0].records = NULL;
              }
            if (streamptr->tsteps[0].recIDs)
              {
                Free(streamptr->tsteps[0].recIDs);
                streamptr->tsteps[0].recIDs = NULL;
              }
          }
        break;
      }
#endif
    default:
      {
        Error("%s support not compiled in (fileID = %d)!", strfiletype(filetype), fileID);
        break;
      }
    }
}

static void
deallocate_sleveltable_t(sleveltable_t *entry)
{
  if (entry->recordID) Free(entry->recordID);
  if (entry->lindex) Free(entry->lindex);
  entry->recordID = NULL;
  entry->lindex = NULL;
}

static void
streamDestroy(stream_t *streamptr)
{
  xassert(streamptr);
  int vlistID = streamptr->vlistID;

  void (*streamCloseDelegate)(stream_t * streamptr, int recordBufIsToBeDeleted)
      = (void (*)(stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_CLOSE_BACKEND).func;

  if (streamptr->filetype != CDI_FILETYPE_UNDEF) streamCloseDelegate(streamptr, 1);

  if (streamptr->record)
    {
      if (streamptr->record->buffer) Free(streamptr->record->buffer);
      Free(streamptr->record);
      streamptr->record = NULL;
    }

  streamptr->filetype = CDI_FILETYPE_UNDEF;
  if (streamptr->filename)
    {
      Free(streamptr->filename);
      streamptr->filename = NULL;
    }

  if (streamptr->vars)
    {
      for (int index = 0; index < streamptr->nvars; index++)
        {
          sleveltable_t *pslev = streamptr->vars[index].recordTable;
          unsigned nsub = streamptr->vars[index].subtypeSize >= 0 ? (unsigned) streamptr->vars[index].subtypeSize : 0U;
          for (size_t isub = 0; isub < nsub; isub++) deallocate_sleveltable_t(pslev + isub);
          if (pslev) Free(pslev);
        }
      Free(streamptr->vars);
      streamptr->vars = NULL;
    }

  if (streamptr->tsteps)
    {
      int maxSteps = streamptr->tstepsNextID;
      for (int index = 0; index < maxSteps; ++index)
        {
          tsteps_t *tstep = &(streamptr->tsteps[index]);
          if (tstep->records) Free(tstep->records);
          if (tstep->recIDs) Free(tstep->recIDs);
          taxisDestroyKernel(&(tstep->taxis));
        }

      Free(streamptr->tsteps);
      streamptr->tsteps = NULL;
    }

#ifdef HAVE_LIBFDB5
  if (streamptr->fdbKeyValueList)
    {
      cdi_fdb_delete_kvlist(streamptr->fdbNumItems, streamptr->fdbKeyValueList);
      streamptr->fdbNumItems = 0;
      streamptr->fdbKeyValueList = NULL;
    }
#endif

  if (vlistID != -1)
    {
      int taxisID = (streamptr->filemode != 'w') ? vlistInqTaxis(vlistID) : -1;
      if (taxisID != -1) taxisDestroy(taxisID);
      void (*mycdiVlistDestroy_)(int, bool) = (void (*)(int, bool)) namespaceSwitchGet(NSSWITCH_VLIST_DESTROY_).func;
      mycdiVlistDestroy_(vlistID, true);
    }

  if (streamptr->jobs) free(streamptr->jobs);
  if (streamptr->jobManager) AsyncWorker_finalize((AsyncManager *) streamptr->jobManager);

  Free(streamptr);
}

static void
streamDestroyP(void *streamptr)
{
  streamDestroy((stream_t *) streamptr);
}

/*
@Function  streamClose
@Title     Close an open dataset

@Prototype  void streamClose(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamClose} closes an open dataset.

@EndFunction
*/
void
streamClose(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  bool lockIO = streamptr->lockIO;
  if (lockIO) CDI_IO_LOCK();

  if (CDI_Debug) Message("streamID = %d filename = %s", streamID, streamptr->filename);
  streamDestroy(streamptr);
  reshRemove(streamID, &streamOps);
  if (CDI_Debug) Message("Removed stream %d from stream list", streamID);

  if (lockIO) CDI_IO_UNLOCK();
}

void
cdiStreamSync_(stream_t *streamptr)
{
  int fileID = streamptr->fileID;
  int filetype = streamptr->filetype;
  int vlistID = streamptr->vlistID;
  int nvars = vlistNvars(vlistID);

  if (fileID == CDI_UNDEFID)
    Warning("File %s not open!", streamptr->filename);
  else if (vlistID == CDI_UNDEFID)
    Warning("Vlist undefined for file %s!", streamptr->filename);
  else if (nvars == 0)
    Warning("No variables defined!");
  else
    {
      if (streamptr->filemode == 'w' || streamptr->filemode == 'a')
        {
          switch (cdiBaseFiletype(filetype))
            {
#ifdef HAVE_LIBNETCDF
            case CDI_FILETYPE_NETCDF:
              {
                void cdf_sync(int ncid);
                if (streamptr->ncmode == 2) cdf_sync(fileID);
                break;
              }
#endif
            default:
              {
                fileFlush(fileID);
                break;
              }
            }
        }
    }
}

/*
@Function  streamSync
@Title     Synchronize an Open Dataset to Disk

@Prototype  void streamSync(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.

@Description
The function @func{streamSync} offers a way to synchronize the disk copy of a dataset with in-memory buffers.

@EndFunction
*/
void
streamSync(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  void (*myStreamSync_)(stream_t * streamptr) = (void (*)(stream_t *)) namespaceSwitchGet(NSSWITCH_STREAM_SYNC).func;
  myStreamSync_(streamptr);
}

int
cdiStreamDefTimestep_(stream_t *streamptr, int tsID)
{
  stream_check_ptr(__func__, streamptr);

  if (CDI_Debug) Message("streamID = %d  tsID = %d", streamptr->self, tsID);

  int vlistID = streamptr->vlistID;
  if (vlistID == CDI_UNDEFID)
    Error("Must not call streamDefTimestep for stream (ID=%d) with (not yet) defined vlist", streamptr->self);

  if (tsID > 0)
    {
      int newtsID = tstepsNewEntry(streamptr);
      if (tsID != newtsID) Error("Internal problem: tsID = %d newtsID = %d", tsID, newtsID);
    }

  int taxisID = vlistInqTaxis(vlistID);
  if (taxisID != CDI_UNDEFID) ptaxisCopy(&streamptr->tsteps[tsID].taxis, taxisPtr(taxisID));

  streamptr->curTsID = tsID;
  streamptr->ntsteps = tsID + 1;

#ifdef HAVE_LIBNETCDF
  int timeIsVarying = vlistHasTime(vlistID);
  if (cdiBaseFiletype(streamptr->filetype) == CDI_FILETYPE_NETCDF && timeIsVarying)
    {
      /* usually points to cdfDefTimestep in serial mode but
       * to cdiPioCdfDefTimestep on servers and to a null-op on
       * clients in client/server mode */
      void (*myCdfDefTimestep)(stream_t * streamptr, int tsID, size_t)
          = (void (*)(stream_t *, int, size_t)) namespaceSwitchGet(NSSWITCH_CDF_DEF_TIMESTEP).func;
      myCdfDefTimestep(streamptr, tsID, 1);
    }
#endif

  cdi_create_records(streamptr, tsID);

  return (int) streamptr->ntsteps;
}

/*
@Function  streamDefTimestep
@Title     Define a timestep

@Prototype int streamDefTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamDefTimestep} defines a timestep of a stream by the identifier tsID.
The identifier tsID is the timestep index starting at 0 for the first timestep.
Before calling this function the functions @func{taxisDefVdate} and @func{taxisDefVtime} should be used
to define the timestamp for this timestep. All calls to write the data refer to this timestep.

@Result
@func{streamDefTimestep} returns the number of expected records of the timestep.

@EndFunction
*/
int
streamDefTimestep(int streamID, int tsID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  if (streamptr->lockIO) CDI_IO_LOCK();

  int (*myStreamDefTimestep_)(stream_t * streamptr, int tsID)
      = (int (*)(stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_DEF_TIMESTEP_).func;
  int status = myStreamDefTimestep_(streamptr, tsID);

  if (streamptr->lockIO) CDI_IO_UNLOCK();

  return status;
}

int
streamInqCurTimestepID(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->curTsID;
}

/*
@Function  streamInqTimestep
@Title     Get timestep information

@Prototype int streamInqTimestep(int streamID, int tsID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.
    @Item  tsID      Timestep identifier.

@Description
The function @func{streamInqTimestep} sets the next timestep to the identifier tsID.
The identifier tsID is the timestep index starting at 0 for the first timestep.
After a call to this function the functions @func{taxisInqVdate} and @func{taxisInqVtime} can be used
to read the timestamp for this timestep. All calls to read the data refer to this timestep.

@Result
@func{streamInqTimestep} returns the number of records of the timestep or 0, if the end of the file is reached.

@EndFunction
*/
int
streamInqTimestep(int streamID, int tsID)
{
  int nrecs = 0;
  stream_t *streamptr = stream_to_pointer(streamID);
  int vlistID = streamptr->vlistID;

  if (tsID < streamptr->ntsteps) streamptr->tsteps[tsID].curRecID = CDI_UNDEFID;  // fix for netCDF
  if (tsID < streamptr->rtsteps)
    {
      streamptr->curTsID = tsID;
      nrecs = streamptr->tsteps[tsID].nrecs;
      streamptr->tsteps[tsID].curRecID = CDI_UNDEFID;
      int taxisID = vlistInqTaxis(vlistID);
      if (taxisID == -1) Error("Timestep undefined for fileID = %d", streamID);
      ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

      return nrecs;
    }

  if (tsID >= streamptr->ntsteps && streamptr->ntsteps != CDI_UNDEFID) return 0;

  int filetype = streamptr->filetype;

  if (CDI_Debug) Message("streamID = %d  tsID = %d  filetype = %d", streamID, tsID, filetype);

  if (streamptr->lockIO) CDI_IO_LOCK();

  switch (cdiBaseFiletype(filetype))
    {
#ifdef HAVE_LIBGRIB
    case CDI_FILETYPE_GRIB:
      {
        switch (streamptr->protocol)
          {
          case CDI_PROTOCOL_FDB: nrecs = fdbInqTimestep(streamptr, tsID); break;

          case CDI_PROTOCOL_ACROSS:  // TODO read from ACROSS
          case CDI_PROTOCOL_OTHER:
          case CDI_PROTOCOL_FILE: nrecs = grbInqTimestep(streamptr, tsID); break;
          }
        break;
      }
#endif
#ifdef HAVE_LIBSERVICE
    case CDI_FILETYPE_SRV:
      {
        nrecs = srvInqTimestep(streamptr, tsID);
        break;
      }
#endif
#ifdef HAVE_LIBEXTRA
    case CDI_FILETYPE_EXT:
      {
        nrecs = extInqTimestep(streamptr, tsID);
        break;
      }
#endif
#ifdef HAVE_LIBIEG
    case CDI_FILETYPE_IEG:
      {
        nrecs = iegInqTimestep(streamptr, tsID);
        break;
      }
#endif
#ifdef HAVE_LIBNETCDF
    case CDI_FILETYPE_NETCDF:
      {
        nrecs = cdfInqTimestep(streamptr, tsID);
        break;
      }
#endif
    default:
      {
        Error("%s support not compiled in!", strfiletype(filetype));
        break;
      }
    }

  if (streamptr->lockIO) CDI_IO_UNLOCK();

  int taxisID = vlistInqTaxis(vlistID);
  if (taxisID == -1) Error("Timestep undefined for fileID = %d", streamID);

  ptaxisCopy(taxisPtr(taxisID), &streamptr->tsteps[tsID].taxis);

  return nrecs;
}

// This function is used in CDO!
SizeType
streamNvals(int streamID)
{
  stream_t *streamptr = stream_to_pointer(streamID);
  return streamptr->numvals;
}

/*
@Function  streamDefVlist
@Title     Define the variable list

@Prototype void streamDefVlist(int streamID, int vlistID)
@Parameter
    @Item  streamID Stream ID, from a previous call to @fref{streamOpenWrite}.
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.

@Description
The function @func{streamDefVlist} defines the variable list of a stream.

To safeguard against errors by modifying the wrong vlist object,
this function makes the passed vlist object immutable.
All further vlist changes have to use the vlist object returned by streamInqVlist().

@EndFunction
*/
void
streamDefVlist(int streamID, int vlistID)
{
  void (*myStreamDefVlist)(int streamID, int vlistID) = (void (*)(int, int)) namespaceSwitchGet(NSSWITCH_STREAM_DEF_VLIST_).func;
  myStreamDefVlist(streamID, vlistID);
}

// The single image implementation of streamDefVlist
void
cdiStreamDefVlist_(int streamID, int vlistID)
{
  stream_t *streamptr = stream_to_pointer(streamID);

  if (streamptr->vlistID == CDI_UNDEFID)
    {
      if (streamptr->lockIO) CDI_IO_LOCK();

      int vlistCopy = vlistDuplicate(vlistID);
      cdiVlistMakeInternal(vlistCopy);
      cdiVlistMakeImmutable(vlistID);
      cdiStreamSetupVlist(streamptr, vlistCopy);

      if (streamptr->lockIO) CDI_IO_UNLOCK();
    }
  else
    Warning("vlist already defined for %s!", streamptr->filename);
}

/*
@Function  streamInqVlist
@Title     Get the variable list

@Prototype int streamInqVlist(int streamID)
@Parameter
    @Item  streamID  Stream ID, from a previous call to @fref{streamOpenRead} or @fref{streamOpenWrite}.

@Description
The function @func{streamInqVlist} returns the variable list of a stream.

@Result
@func{streamInqVlist} returns an identifier to the variable list.

@EndFunction
*/
int
streamInqVlist(int streamID)
{
  stream_t *s = stream_to_pointer(streamID);
  return s->vlistID;
}

void
streamDefShuffle(int streamID, int shuffle)
{
  stream_t *s = stream_to_pointer(streamID);
  if (s->shuffle != shuffle)
    {
      s->shuffle = shuffle;
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}

void
streamDefFilter(int streamID, int filterId, int numParams, const int *params)
{
  stream_t *s = stream_to_pointer(streamID);
  if ((int) s->filterId != filterId)
    {
      if (numParams > (int) s->maxParams) Error("Too many filter parameter %d (max=%zu)!", numParams, s->maxParams);
      s->filterId = filterId;
      s->numParams = numParams;
      for (int i = 0; i < numParams; ++i) s->params[i] = params[i];
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}

void
streamDefCompType(int streamID, int comptype)
{
  stream_t *s = stream_to_pointer(streamID);
  if (s->comptype != comptype)
    {
      s->comptype = comptype;
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}

void
streamDefCompLevel(int streamID, int complevel)
{
  stream_t *s = stream_to_pointer(streamID);
  if (s->complevel != complevel)
    {
      s->complevel = complevel;
      reshSetStatus(streamID, &streamOps, RESH_DESYNC_IN_USE);
    }
}

int
streamInqCompType(int streamID)
{
  stream_t *s = stream_to_pointer(streamID);
  return s->comptype;
}

int
streamInqCompLevel(int streamID)
{
  stream_t *s = stream_to_pointer(streamID);
  return s->complevel;
}

int
streamInqFileID(int streamID)
{
  stream_t *s = (stream_t *) reshGetVal(streamID, &streamOps);
  return s->fileID;
}

void
cdiDefAccesstype(int streamID, int type)
{
  stream_t *s = (stream_t *) reshGetVal(streamID, &streamOps);

  if (s->accesstype == CDI_UNDEFID)
    {
      s->accesstype = type;
    }
  else if (s->accesstype != type)
    Error("Changing access type from %s not allowed!", s->accesstype == TYPE_REC ? "REC to VAR" : "VAR to REC");
}

int
cdiInqAccesstype(int streamID)
{
  stream_t *s = (stream_t *) reshGetVal(streamID, &streamOps);
  return s->accesstype;
}

static int
streamTxCode(void *s)
{
  (void) s;
  return STREAM;
}

void
cdiStreamSetupVlist(stream_t *s, int vlistID)
{
  void (*myStreamSetupVlist)(stream_t * s, int vlistID)
      = (void (*)(stream_t *, int)) namespaceSwitchGet(NSSWITCH_STREAM_SETUP_VLIST).func;
  myStreamSetupVlist(s, vlistID);
}

void
cdiStreamSetupVlist_(stream_t *streamptr, int vlistID)
{
  streamptr->vlistID = vlistID;
  int nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
    {
      int gridID = vlistInqVarGrid(vlistID, varID);
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int tilesetID = vlistInqVarSubtype(vlistID, varID);
      stream_new_var(streamptr, gridID, zaxisID, tilesetID);
      if (streamptr->have_missval) vlistDefVarMissval(vlistID, varID, vlistInqVarMissval(vlistID, varID));
    }

  if (streamptr->filemode == 'w')
    {
      tstepsNewEntry(streamptr);  // timestep 0
      int vlistIDw = streamptr->vlistID;
      int timeIsVarying = vlistHasTime(vlistIDw);
      if (timeIsVarying)
        {
          int taxisID = vlistInqTaxis(vlistIDw);
          if (taxisID == CDI_UNDEFID)
            {
              Warning("taxisID undefined for fileID = %d! Using absolute time axis.", streamptr->self);
              taxisID = taxisCreate(TAXIS_ABSOLUTE);
              vlistDefTaxis(vlistIDw, taxisID);
            }

#ifdef HAVE_LIBNETCDF
          if (taxisInqType(taxisID) == TAXIS_RELATIVE)
            if (cdiBaseFiletype(streamptr->filetype) == CDI_FILETYPE_NETCDF)
              {
                const taxis_t *taxisptr = taxisPtr(taxisID);
                if (cdiDateTime_isNull(taxisptr->rDateTime))
                  {
                    int vdate = taxisInqVdate(taxisID);
                    if (vdate == 0) vdate = 10101;
                    taxisDefRdate(taxisID, vdate);
                  }
              }
#endif
          ptaxisCopy(&streamptr->tsteps[0].taxis, taxisPtr(taxisID));
        }

      switch (cdiBaseFiletype(streamptr->filetype))
        {
#ifdef HAVE_LIBNETCDF
        case CDI_FILETYPE_NETCDF:
          {
            /* calls cdfDefCoordinateVars in serial mode but
             * cdiPioClientStreamNOP (i.e. nothing) on client ranks
             * and cdiPioServerCdfDefVars on server ranks in parallel mode*/
            void (*myCdfDefVars)(stream_t * streamptr) = (void (*)(stream_t *)) namespaceSwitchGet(NSSWITCH_CDF_STREAM_SETUP).func;
            myCdfDefVars(streamptr);
          }
          break;
#endif
#ifdef HAVE_LIBGRIB
        case CDI_FILETYPE_GRIB: gribContainersNew(streamptr); break;
#endif
        default:;
        }
    }
}

void
cdiStreamGetIndexList(unsigned numIDs, int *IDs)
{
  reshGetResHListOfType(numIDs, IDs, &streamOps);
}

int
streamInqNvars(int streamID)
{
  stream_t *s = (stream_t *) reshGetVal(streamID, &streamOps);
  return s->nvars;
}

static int
streamCompareP(void *streamptr1, void *streamptr2)
{
  stream_t *s1 = (stream_t *) streamptr1;
  stream_t *s2 = (stream_t *) streamptr2;
  enum
  {
    differ = -1,
    equal = 0,
  };

  xassert(s1);
  xassert(s2);

  if (s1->filetype != s2->filetype) return differ;
  if (s1->byteorder != s2->byteorder) return differ;
  if (s1->comptype != s2->comptype) return differ;
  if (s1->complevel != s2->complevel) return differ;

  if (s1->filename)
    {
      if (!str_is_equal(s1->filename, s2->filename)) return differ;
    }
  else if (s2->filename)
    return differ;

  return equal;
}

void
streamPrintP(void *streamptr, FILE *fp)
{
  stream_t *sp = (stream_t *) streamptr;

  if (!sp) return;

  fprintf(fp,
          "#\n"
          "# streamID %d\n"
          "#\n"
          "self          = %d\n"
          "accesstype    = %d\n"
          "accessmode    = %d\n"
          "filetype      = %d\n"
          "byteorder     = %d\n"
          "fileID        = %d\n"
          "filemode      = %d\n"
          "filename      = %s\n"
          "nrecs         = %d\n"
          "nvars         = %d\n"
          "varsAllocated = %d\n"
          "curTsID       = %d\n"
          "rtsteps       = %d\n"
          "ntsteps       = %ld\n"
          "tstepsTableSize= %d\n"
          "tstepsNextID  = %d\n"
          "ncmode        = %d\n"
          "vlistID       = %d\n"
          "globalatts    = %d\n"
          "localatts     = %d\n"
          "unreduced     = %d\n"
          "sortname      = %d\n"
          "have_missval  = %d\n"
          "ztype         = %d\n"
          "zlevel        = %d\n",
          sp->self, sp->self, sp->accesstype, sp->accessmode, sp->filetype, sp->byteorder, sp->fileID, sp->filemode, sp->filename,
          sp->nrecs, sp->nvars, sp->varsAllocated, sp->curTsID, sp->rtsteps, sp->ntsteps, sp->tstepsTableSize, sp->tstepsNextID,
          sp->ncmode, sp->vlistID, sp->globalatts, sp->localatts, sp->unreduced, sp->sortname, sp->have_missval, sp->comptype,
          sp->complevel);
}

enum
{
  streamNint = 10,
};

static int
streamGetPackSize(void *voidP, void *context)
{
  stream_t *streamP = (stream_t *) voidP;
  int packBufferSize = serializeGetSize(streamNint, CDI_DATATYPE_INT, context) + serializeGetSize(2, CDI_DATATYPE_UINT32, context)
                       + serializeGetSize((int) strlen(streamP->filename) + 1, CDI_DATATYPE_TXT, context)
                       + serializeGetSize(1, CDI_DATATYPE_FLT64, context);
  return packBufferSize;
}

static void
streamPack(void *streamptr, void *packBuffer, int packBufferSize, int *packBufferPos, void *context)
{
  stream_t *streamP = (stream_t *) streamptr;
  int intBuffer[streamNint];

  intBuffer[0] = streamP->self;
  intBuffer[1] = streamP->filetype;
  intBuffer[2] = (int) strlen(streamP->filename) + 1;
  intBuffer[3] = streamP->vlistID;
  intBuffer[4] = streamP->byteorder;
  intBuffer[5] = streamP->comptype;
  intBuffer[6] = streamP->complevel;
  intBuffer[7] = streamP->unreduced;
  intBuffer[8] = streamP->sortname;
  intBuffer[9] = streamP->have_missval;

  serializePack(intBuffer, streamNint, CDI_DATATYPE_INT, packBuffer, packBufferSize, packBufferPos, context);
  uint32_t d = cdiCheckSum(CDI_DATATYPE_INT, streamNint, intBuffer);
  serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);

  serializePack(&CDI_Default_Missval, 1, CDI_DATATYPE_FLT64, packBuffer, packBufferSize, packBufferPos, context);
  serializePack(streamP->filename, intBuffer[2], CDI_DATATYPE_TXT, packBuffer, packBufferSize, packBufferPos, context);
  d = cdiCheckSum(CDI_DATATYPE_TXT, intBuffer[2], streamP->filename);
  serializePack(&d, 1, CDI_DATATYPE_UINT32, packBuffer, packBufferSize, packBufferPos, context);
}

struct streamAssoc
streamUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context)
{
  int intBuffer[streamNint];
  uint32_t d;
  char filename[CDI_MAX_NAME];

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, intBuffer, streamNint, CDI_DATATYPE_INT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
  xassert(cdiCheckSum(CDI_DATATYPE_INT, streamNint, intBuffer) == d);

  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &CDI_Default_Missval, 1, CDI_DATATYPE_FLT64, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &filename, intBuffer[2], CDI_DATATYPE_TXT, context);
  serializeUnpack(unpackBuffer, unpackBufferSize, unpackBufferPos, &d, 1, CDI_DATATYPE_UINT32, context);
  xassert(d == cdiCheckSum(CDI_DATATYPE_TXT, intBuffer[2], filename));
  int targetStreamID = namespaceAdaptKey(intBuffer[0], originNamespace),
      streamID = streamOpenID(filename, 'w', intBuffer[1], targetStreamID);
  xassert(streamID >= 0 && targetStreamID == streamID);
  streamDefByteorder(streamID, intBuffer[4]);
  streamDefCompType(streamID, intBuffer[5]);
  streamDefCompLevel(streamID, intBuffer[6]);
  stream_t *streamptr = stream_to_pointer(streamID);
  streamptr->unreduced = intBuffer[7];
  streamptr->sortname = intBuffer[8];
  streamptr->have_missval = intBuffer[9];
  struct streamAssoc retval = { streamID, intBuffer[3] };
  return retval;
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
