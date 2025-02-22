#ifndef _STREAM_EXT_H
#define _STREAM_EXT_H

#ifndef _EXTRA_H
#include "extra.h"
#endif

int extInqContents(stream_t *streamptr);
int extInqTimestep(stream_t *streamptr, int tsID);

int extInqRecord(stream_t *streamptr, int *varID, int *levelID);
void extDefRecord(stream_t *streamptr);
void extCopyRecord(stream_t *streamptr2, stream_t *streamptr1);
void ext_read_record(stream_t *streamptr, int memtype, void *data, size_t *numMissVals);
void ext_write_record(stream_t *streamptr, int memtype, const void *data);

void extReadVarDP(stream_t *streamptr, int varID, double *data, size_t *numMissVals);
void extWriteVarDP(stream_t *streamptr, int varID, const double *data);

void extReadVarSliceDP(stream_t *streamptr, int varID, int levelID, double *data, size_t *numMissVals);
void extWriteVarSliceDP(stream_t *streamptr, int varID, int levelID, const double *data);

#endif /* _STREAM_EXT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
