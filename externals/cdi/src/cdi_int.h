#ifndef CDI_INT_H
#define CDI_INT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>

#include "cdi.h"
#include "cdf_config.h"

// Base file types

#define CDI_FILETYPE_GRIB 100    // File type GRIB
#define CDI_FILETYPE_NETCDF 101  // File type NetCDF

// dummy use of unused parameters to silence compiler warnings
#ifndef UNUSED
#define UNUSED(x) (void) x
#endif

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = (char *) Malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

char *str_to_lower(char *str);
bool strStartsWith(const char *vstr, const char *cstr);

static inline bool
str_is_equal(const char *x, const char *y)
{
  return (*x == *y) && strcmp(x, y) == 0;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

#ifndef ERROR_H
#include "error.h"
#endif
#ifndef _BASETIME_H
#include "basetime.h"
#endif
#ifndef JULIAN_DATE_H
#include "julian_date.h"
#endif
#ifndef TAXIS_H
#include "taxis.h"
#endif
#ifndef CDI_LIMITS_H
#include "cdi_limits.h"
#endif
#ifndef _SERVICE_H
#include "service.h"
#endif
#ifndef _EXTRA_H
#include "extra.h"
#endif
#ifndef _IEG_H
#include "ieg.h"
#endif
#ifndef RESOURCE_HANDLE_H
#include "resource_handle.h"
#endif

#define check_parg(arg) \
  if (arg == 0) Warning("Argument '" #arg "' not allocated!")

#ifdef __xlC__ /* performance problems on IBM */
#ifndef DBL_IS_NAN
#define DBL_IS_NAN(x) ((x) != (x))
#endif
#else
#ifndef DBL_IS_NAN
#if defined(HAVE_DECL_ISNAN)
#define DBL_IS_NAN(x) (isnan(x))
#elif defined(FP_NAN)
#define DBL_IS_NAN(x) (fpclassify(x) == FP_NAN)
#else
#define DBL_IS_NAN(x) ((x) != (x))
#endif
#endif
#endif

#ifndef DBL_IS_EQUAL
//#define DBL_IS_EQUAL(x,y) (!(x < y || y < x))
#define DBL_IS_EQUAL(x, y) (DBL_IS_NAN(x) || DBL_IS_NAN(y) ? (DBL_IS_NAN(x) && DBL_IS_NAN(y)) : !(x < y || y < x))
#endif

#ifndef IS_EQUAL
#define IS_NOT_EQUAL(x, y) (x < y || y < x)
#define IS_EQUAL(x, y) (!IS_NOT_EQUAL(x, y))
#endif

enum
{
  TYPE_REC,
  TYPE_VAR,
};

enum
{
  MEMTYPE_DOUBLE = 1,
  MEMTYPE_FLOAT,
};

typedef struct
{
  void *buffer;       // gribapi, cgribex
  size_t buffersize;  // gribapi, cgribex
  off_t position;     // file position
  int param;
  int ilevel;
  int vdate;
  int vtime;
  int gridID;
  int varID;
  int levelID;
  int prec;       // ext, srv
  void *objectp;  // pointer to ieg, ext, srv or cgribex objects
} Record;

// data structure specifying tile-related meta-data. structure contains "-1" if this is no tile-variable.
typedef struct
{
  int tileindex, totalno_of_tileattr_pairs, tileClassification, numberOfTiles, numberOfAttributes, attribute;
} var_tile_t;

typedef struct
{
  short perturbationNumber;
  short typeOfGeneratingProcess;
} VarScanKeys;

static inline void
varScanKeysInit(VarScanKeys *s)
{
  memset(s, 0, sizeof(VarScanKeys));
}

static inline bool
varScanKeysIsEqual(const VarScanKeys *s1, const VarScanKeys *s2)
{
  return memcmp(s1, s2, sizeof(VarScanKeys)) == 0;
}

typedef struct
{
  off_t position;
  size_t size;
  size_t gridsize;
  int zip;
  int param;
  int ilevel;
  int ilevel2;
  int ltype;
  short tsteptype;
  short varID;
  int levelID;
  short used;
  char varname[32];  // needed for grib decoding with GRIB_API
  VarScanKeys scanKeys;
  var_tile_t tiles;  // tile-related meta-data, currently for GRIB-API only.
#ifdef HAVE_LIBFDB5
  void *fdbItem;
#endif
} record_t;

typedef struct
{
  int *recIDs;  // IDs of non constant records
  record_t *records;
  int recordSize;   // number of allocated records
  int nrecs;        // number of used records
                    // tsID=0 nallrecs
                    // tsID>0 number of non constant records
  int nallrecs;     // number of all records
  int curRecID;     // current record ID
  int ncStepIndex;  // NetCDF timestep index
  off_t position;   // timestep file position
  taxis_t taxis;
  bool next;
} tsteps_t;

typedef struct
{
  int nlevs;
  int subtypeIndex;  // corresponding tile in subtype_t structure (subtype->self)
  int *recordID;     // record IDs: [nlevs]
  int *lindex;       // level index
} sleveltable_t;

typedef struct
{
  sleveltable_t *recordTable;  // record IDs for each subtype
  int ncvarid;
  int subtypeSize;
  bool defmiss;  // true: if missval is defined in file
  bool isUsed;

  int gridID;
  int zaxisID;
  int tsteptype;  // TSTEP_*
  int subtypeID;  // subtype ID, e.g. for tile-related meta-data (currently for GRIB-API only).
} svarinfo_t;

typedef struct
{
  int ilev;
  int mlev;
  int ilevID;
  int mlevID;
} VCT;

#ifdef HAVE_LIBNETCDF
enum cdfIDIdx
{
  CDF_DIMID_E,  // 3rd dimID of cube sphere grid (len=6)
  CDF_DIMID_X,
  CDF_DIMID_Y,
  CDF_DIMID_RP,  // reducedPoints
  CDF_VARID_X,
  CDF_VARID_Y,
  CDF_VARID_RP,  // reducedPoints
  CDF_VARID_A,
  CDF_SIZE_ncIDs,
};

typedef struct
{
  int ncIDs[CDF_SIZE_ncIDs];
  int gridID;
  long start;
  long count;
} ncgrid_t;
#endif

typedef struct
{
  int self;
  int accesstype;  // TYPE_REC or TYPE_VAR
  int accessmode;
  int filetype;
  int byteorder;
  int fileID;
  int filemode;
  int nrecs;  // number of records
  SizeType numvals;
  char *filename;
  Record *record;
  CdiQuery *query;
  svarinfo_t *vars;
  int nvars;  // number of variables
  int varsAllocated;
  int curTsID;   // current timestep ID
  int rtsteps;   // number of tsteps accessed
  long ntsteps;  // number of tsteps : only set if all records accessed
  int maxSteps;  // max. number of timesteps (needed for CDI_FILETYPE_NCZARR)
  tsteps_t *tsteps;
  int tstepsTableSize;
  int tstepsNextID;
  basetime_t basetime;
  int ncmode;
  int vlistID;
#ifdef HAVE_LIBNETCDF
  int nc_complex_float_id;
  int nc_complex_double_id;
  ncgrid_t ncgrid[MAX_GRIDS_PS];
  int zaxisID[MAX_ZAXES_PS];  // Warning: synchronous array to vlist_to_pointer(vlistID)->zaxisIDs
  int nczvarID[MAX_ZAXES_PS];
  int ncNumDims;
  int ncDimID[MAX_DIMS_PS];
  size_t ncDimLen[MAX_DIMS_PS];
  VCT vct;
  size_t chunkSizeTdim;
  size_t chunkSizeZdim;
#endif
  long maxGlobalRecs;
  int globalatts;
  int localatts;
  int unreduced;
  int have_missval;
  // netcdf4/HDF5 filter
  unsigned int filterId;
  size_t numParams;
  size_t maxParams;
  unsigned int params[8];

  int comptype;   // compression type
  int complevel;  // compression level
  bool sortname;
  bool sortparam;
  void *gribContainers;

  int numWorker;
  int nextGlobalRecId;
  int cachedTsID;
  void *jobs;
  void *jobManager;

  int protocol;
  void *protocolData;
} stream_t;

// Length of optional keyword/value pair list
#define MAX_OPT_GRIB_ENTRIES 500

enum cdi_convention
{
  CDI_CONVENTION_ECHAM,
  CDI_CONVENTION_CF
};

// Data type specification for optional key/value pairs (GRIB)
typedef enum
{
  t_double = 0,
  t_int = 1
} key_val_pair_datatype;

// Data structure holding optional key/value pairs for GRIB
typedef struct
{
  char *keyword;  // keyword string
  bool update;
  key_val_pair_datatype data_type;  // data type of this key/value pair
  double dbl_val;                   // double value (data_type == t_double)
  int int_val;                      // integer value (data_type == t_int)
  int subtype_index;                // tile index for this key-value pair
} opt_key_val_pair_t;

// enum for differenciating between the different times that we handle
typedef enum
{
  kCdiTimeType_referenceTime,
  kCdiTimeType_startTime,
  kCdiTimeType_endTime
} CdiTimeType;

#define CDI_FILETYPE_UNDEF -1  // Unknown/not yet defined file type

extern int cdiDebugExt;
extern int CDI_Debug;  // If set to 1, debuggig (default 0)
extern int CDI_Recopt;
extern bool CDI_gribapi_debug;
extern bool CDI_gribapi_grib1;
extern double CDI_Default_Missval;
extern double CDI_Grid_Missval;
extern int CDI_Default_InstID;
extern int CDI_Default_ModelID;
extern int CDI_Default_TableID;
extern int cdiDefaultLeveltype;
extern int CDI_Default_Calendar;
// extern int cdiNcMissingValue;
extern int CDI_Netcdf_Chunksizehint;
extern int CDI_ChunkType;
extern int CDI_Test;
extern int CDI_Split_Ltype105;
extern int cdiDataUnreduced;
extern int cdiSortName;
extern int cdiSortParam;
extern int cdiHaveMissval;
extern bool CDI_Ignore_Att_Coordinates;
extern bool CDI_Coordinates_Lon_Lat;
extern bool CDI_Ignore_Valid_Range;
extern int CDI_Skip_Records;
extern const char *CDI_GRIB1_Template;
extern const char *CDI_GRIB2_Template;
extern int CDI_Convention;
extern int CDI_Inventory_Mode;
extern int CDI_Query_Abort;
extern int CDI_Version_Info;
extern int CDI_Convert_Cubesphere;
extern int CDI_Read_Cell_Corners;
extern int CDI_CMOR_Mode;
extern int CDI_Reduce_Dim;
extern int CDI_Shuffle;
extern size_t CDI_Netcdf_Hdr_Pad;
extern size_t CDI_Chunk_Cache;
extern size_t CDI_Chunk_Cache_Max;
extern bool CDI_Netcdf_Lazy_Grid_Load;
extern int STREAM_Debug;

extern char *cdiPartabPath;
extern int cdiPartabIntern;
extern const resOps streamOps;

static inline stream_t *
stream_to_pointer(int idx)
{
  return (stream_t *) reshGetVal(idx, &streamOps);
}

static inline void
stream_check_ptr(const char *caller, stream_t *streamptr)
{
  if (streamptr == NULL) Errorc("stream undefined!");
}

int streamInqFileID(int streamID);

void gridDefHasDims(int gridID, int hasdims);
int gridInqHasDims(int gridID);
int zaxisInqLevelID(int zaxisID, double level);

void streamCheckID(const char *caller, int streamID);

void streamDefineTaxis(int streamID);

int streamsNewEntry(int filetype);
void streamsInitEntry(int streamID);
void cdiStreamSetupVlist(stream_t *streamptr, int vlistID);
// default implementation of the overridable function
void cdiStreamSetupVlist_(stream_t *streamptr, int vlistID);
int stream_new_var(stream_t *streamptr, int gridID, int zaxisID, int tilesetID);

int tstepsNewEntry(stream_t *streamptr);

const char *strfiletype(int filetype);

void cdi_generate_vars(stream_t *streamptr);

void vlist_check_contents(int vlistID);

void cdi_create_records(stream_t *streamptr, int tsID);

void streamFCopyRecord(stream_t *streamptr2, stream_t *streamptr1, const char *container_name);

int recordNewEntry(stream_t *streamptr, int tsID);

void cdi_create_timesteps(int numTimesteps, stream_t *streamptr);

void recordInitEntry(record_t *record);

void cdiCheckZaxis(int zaxisID);

void cdiDefAccesstype(int streamID, int type);
int cdiInqAccesstype(int streamID);

int getByteswap(int byteorder);

void cdiStreamGetIndexList(unsigned numIDs, int IDs[]);

void cdiInitialize(void);

char *cdiEscapeSpaces(const char *string);
char *cdiUnescapeSpaces(const char *string, const char **outStringEnd);

enum
{
  CDI_UNIT_PA = 1,
  CDI_UNIT_HPA,
  CDI_UNIT_MM,
  CDI_UNIT_CM,
  CDI_UNIT_DM,
  CDI_UNIT_M,
};

struct streamAssoc
{
  int streamID, vlistID;
};

struct streamAssoc streamUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context);

int cdiStreamOpenDefaultDelegate(const char *filename, char filemode, int filetype, stream_t *streamptr,
                                 int recordBufIsToBeCreated);

int streamOpenID(const char *filename, char filemode, int filetype, int resH);

void cdiStreamDefVlist_(int streamID, int vlistID);

int cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data, SizeType nmiss);

void cdiStreamWriteVarChunk_(int streamID, int varID, int memtype, const int rect[][2], const void *data, SizeType nmiss);
void cdiStreamCloseDefaultDelegate(stream_t *streamptr, int recordBufIsToBeDeleted);

int cdiStreamDefTimestep_(stream_t *streamptr, int tsID);

void cdiStreamSync_(stream_t *streamptr);

const char *cdiUnitNamePtr(int cdi_unit);

enum
{
  // 8192 is known to work on most systems (4096 isn't on Alpha)
  commonPageSize = 8192,
};

size_t cdiGetPageSize(bool largePageAlign);

void zaxisGetIndexList(int nzaxis, int *zaxisIndexList);

// clang-format off

#ifdef __cplusplus
extern "C" {
#endif

// functions used in CDO !!!

void cdiDefTableID(int tableID);

void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);

static inline
void cdi_check_gridsize_int_limit(const char *format, SizeType gridsize)
{
  if (gridsize > INT_MAX) Error("%s format grid size (%zu) limit exceeded (%zu)!", format, gridsize, INT_MAX);
}

bool cdiFiletypeIsExse(int filetype);
int cdiBaseFiletype(int filetype);

#ifdef __cplusplus
}
#endif

// clang-format on

#endif /* CDI_INT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
