#ifndef _STREAM_CDF_H
#define _STREAM_CDF_H

#include "cdi_int.h"
#include "taxis.h"
#include "grid.h"

enum
{
  POSITIVE_UP = 1,
  POSITIVE_DOWN = 2,
};

enum
{
  CDF_MAX_TIME_UNIT_STR     /* maximum length of time unit string */
  = TAXIS_MAX_UNIT_STR_LEN  /* longest result from tunitNamePtr */
    + 7                     /* room for " since " */
    + 7 + 1 + 2 + 1 + 2     /* room for year with 7 digits,
                             * dashes and 2 digits for month and day */
    + 1 + 2 + 1 + 2 + 1 + 2 /* room for " " and 2 digit hour, minute,
                             *                          second */
    + 1                     /* and terminating '\0' */
};

int cdfDefVar(stream_t *streamptr, int varID);
void cdfDefCoordinateVars(stream_t *streamptr);
void cdfDefTimestep(stream_t *streamptr, int tsID, size_t valCount);
int cdfInqTimestep(stream_t *streamptr, int tsID);
int cdfInqContents(stream_t *streamptr);

void cdfEndDef(stream_t *streamptr);
void cdfDefRecord(stream_t *streamptr);

void cdfCopyRecord(stream_t *streamptr2, stream_t *streamptr1);

void cdfDefineAttributes(int filetype, int vlistID, int varID, int fileID, int ncvarID);

void cdf_read_record(stream_t *streamptr, int memtype, void *data, size_t *numMissVals);
void cdf_write_record(stream_t *streamptr, int memtype, const void *data, size_t numMissVals);

void cdf_read_var(stream_t *streamptr, int varID, int memtype, void *data, size_t *numMissVals);
void cdf_write_var(stream_t *streamptr, int varID, int memtype, const void *data, size_t numMissVals);

void cdf_read_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, void *data, size_t *numMissVals);
void cdf_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, size_t numMissVals);

void cdf_write_var_chunk(stream_t *streamptr, int varID, int memtype, const int rect[][2], const void *data, size_t numMissVals);

void cdfDefVarDeflate(int ncid, int ncvarid, int shuffle, int deflateLevel);
void cdfDefTime(stream_t *streamptr);

void cdf_scale_add(size_t size, double *data, double addoffset, double scalefactor);

int cdfDefDatatype(int datatype, stream_t *streamptr);

void cdf_create_records(stream_t *streamptr, int tsID);

#define ChunkSizeMax 65536
#define ChunkSizeLim 16777216
size_t calc_chunksize_x(int chunkType, int chunkSize, size_t xsize, bool yIsUndefined);
size_t calc_chunksize_y(int chunkType, size_t gridsize, size_t xsize, size_t ysize);

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
