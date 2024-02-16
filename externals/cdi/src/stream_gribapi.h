#ifndef STREAM_GRIBAPI_H
#define STREAM_GRIBAPI_H

#ifdef HAVE_LIBGRIB_API

#include "gribapi.h"

int fdbScanTimesteps(stream_t *streamptr);

int gribapiScanTimestep1(stream_t *streamptr);
int gribapiScanTimestep2(stream_t *streamptr);
int gribapiScanTimestep(stream_t *streamptr);

int gribapiDecode(int memType, void *gribbuffer, size_t gribsize, void *data, size_t datasize, int unreduced, size_t *nmiss,
                  double missval);

size_t gribapiEncode(int memType, int varID, int levelID, int vlistID, int gridID, int zaxisID, CdiDateTime vDateTime,
                     int tsteptype, int numavg, size_t datasize, const void *data, size_t nmiss, void **gribbuffer,
                     size_t *gribbuffersize, int ljpeg, void *gribContainer);

int gribapiGetScanningMode(grib_handle *gh);
void gribapiSetScanningMode(grib_handle *gh, int scanningMode);

void gribapiChangeParameterIdentification(grib_handle *gh, int code, int ltype, int lev);

#endif

#endif /* STREAM_GRIBAPI_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
