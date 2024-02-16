#ifndef TAXIS_H
#define TAXIS_H

#include <stdbool.h>
#include "cdi.h"

#ifndef RESOURCE_HANDLE_H
#include "resource_handle.h"
#endif

typedef struct
{
  int self;
  int datatype;  // datatype
  int type;      // time type
  int calendar;
  int unit;  // time units
  int numavg;
  CdiDateTime sDateTime;     // start date/time
  CdiDateTime vDateTime;     // verification date/time
  CdiDateTime rDateTime;     // reference date/time
  CdiDateTime fDateTime;     // forecast reference date/time
  CdiDateTime vDateTime_lb;  // lower bounds of verification date/time
  CdiDateTime vDateTime_ub;  // upper bounds of verification date/time
  double fc_period;          // forecast time period
  int fc_unit;               // forecast time unit
  char *name;
  char *longname;
  char *units;
  bool climatology;
  bool hasBounds;
} taxis_t;

//      taxisInqSdatetime: Get the start date/time
CdiDateTime taxisInqSdatetime(int taxisID);

void ptaxisInit(taxis_t *taxis);
void ptaxisCopy(taxis_t *dest, taxis_t *source);
taxis_t *taxisPtr(int taxisID);
void cdi_set_forecast_period(double timevalue, taxis_t *taxis);
CdiDateTime cdi_decode_timeval(double timevalue, const taxis_t *taxis);
double cdi_encode_timeval(CdiDateTime datetime, taxis_t *taxis);

void ptaxisDefDatatype(taxis_t *taxisptr, int datatype);
void ptaxisDefName(taxis_t *taxisptr, const char *name);
void ptaxisDefLongname(taxis_t *taxisptr, const char *longname);
void ptaxisDefUnits(taxis_t *taxisptr, const char *units);
char *ptaxisAllocUnits(taxis_t *taxisptr, size_t len);
void taxisDestroyKernel(taxis_t *taxisptr);
#ifndef SX
extern const resOps taxisOps;
#endif

int taxisUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                int checkForSameID);

enum
{
  TAXIS_MAX_UNIT_STR_LEN = 9
};

#endif /* TAXIS_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
