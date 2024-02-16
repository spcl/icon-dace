#ifndef STREAM_CGRIBEX_H
#define STREAM_CGRIBEX_H

void *cgribexNew(void);
void cgribexDelete(void *cgribexp);

int cgribexScanTimestep1(stream_t *streamptr);
int cgribexScanTimestep2(stream_t *streamptr);
int cgribexScanTimestep(stream_t *streamptr);

int cgribexDecode(int memtype, void *cgribexp, void *gribbuffer, size_t gribsize, void *data, size_t datasize, int unreduced,
                  size_t *nmiss, double missval);

size_t cgribexEncode(int memtype, int varID, int levelID, int vlistID, int gridID, int zaxisID, CdiDateTime vDateTime,
                     int tsteptype, int numavg, size_t datasize, const void *data, size_t nmiss, void *gribbuffer,
                     size_t gribbuffersize);

void *cgribex_handle_new_from_meassage(void *gribbuffer, size_t recsize);
void cgribex_handle_delete(void *gh);

void cgribexChangeParameterIdentification(void *gh, int code, int ltype, int lev);

#endif /* STREAM_CGRIBEX_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
