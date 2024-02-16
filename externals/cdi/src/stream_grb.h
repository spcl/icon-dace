#ifndef STREAM_GRB_H
#define STREAM_GRB_H

double zaxis_units_to_centimeter(int zaxisID);
double zaxis_units_to_meter(int zaxisID);
bool zaxis_units_is_Pa(int zaxisID);

void ensureBufferSize(size_t requiredSize, size_t *curSize, void **buffer);
int grbDecompress(size_t recsize, size_t *buffersize, void **gribbuffer);

static inline bool
gribbyte_get_bit(int number, int bit)
{
  return (bool) ((number >> (8 - bit)) & 1);
}
static inline void
gribbyte_set_bit(int *number, int bit)
{
  *number |= 1 << (8 - bit);
}
static inline void
gribbyte_clear_bit(int *number, int bit)
{
  *number &= ~(1 << (8 - bit));
}

int grbBitsPerValue(int datatype);

int fdbInqContents(stream_t *streamptr);
int grbInqContents(stream_t *streamptr);
int fdbInqTimestep(stream_t *streamptr, int tsID);
int grbInqTimestep(stream_t *streamptr, int tsID);

int grbInqRecord(stream_t *streamptr, int *varID, int *levelID);
void grbDefRecord(stream_t *streamptr);
void grb_read_record(stream_t *streamptr, int memtype, void *data, size_t *nmiss);
void grb_write_record(stream_t *streamptr, int memtype, const void *data, size_t nmiss);
void grbCopyRecord(stream_t *streamptr2, stream_t *streamptr1);

void grb_read_var(stream_t *streamptr, int varID, int memtype, void *data, size_t *nmiss);
void grb_write_var(stream_t *streamptr, int varID, int memtype, const void *data, size_t nmiss);

void grb_read_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, void *data, size_t *nmiss);
void grb_write_var_slice(stream_t *streamptr, int varID, int levelID, int memtype, const void *data, size_t nmiss);

int grib1ltypeToZaxisType(int grib_ltype);
int grib2ltypeToZaxisType(int grib_ltype);

int zaxisTypeToGrib1ltype(int zaxistype);
int zaxisTypeToGrib2ltype(int zaxistype);

int grbGetGridtype(int *gridID, size_t gridsize, bool *gridIsRotated, bool *gridIsCurvilinear);

struct cdiGribParamChange
{
  int code, ltype, lev;
  bool active;
};

struct cdiGribScanModeChange
{
  int value;
  bool active;
};

extern struct cdiGribParamChange cdiGribChangeParameterID;
extern struct cdiGribScanModeChange cdiGribDataScanningMode;

// Used in CDO
void streamGrbChangeParameterIdentification(int code, int ltype, int lev);
void streamGrbDefDataScanningMode(int scanmode);

#endif /* STREAM_GRB_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
