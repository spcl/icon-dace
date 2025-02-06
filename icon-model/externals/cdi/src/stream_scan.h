#ifndef STREAM_SCAN_H
#define STREAM_SCAN_H

#include "cdi_int.h"

void streamScanResizeRecords1(stream_t *streamptr);
int streamScanInitRecords2(stream_t *streamptr);
int streamScanInitRecords(stream_t *streamptr, int tsID);
void streamScanTimeConstAdjust(stream_t *streamptr, const taxis_t *taxis);
void streamScanTsFixNtsteps(stream_t *streamptr, off_t recpos);

#endif /* STREAM_SCAN_H */
