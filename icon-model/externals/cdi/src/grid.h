#ifndef GRID_H
#define GRID_H

#include "cdi.h"
#include <stdbool.h>

#include "cdi_key.h"
#include "cdi_att.h"
#include "resource_handle.h"

typedef unsigned char mask_t;

typedef struct grid_t grid_t;

enum gridPropInq
{
  GRID_PROP_MASK,
  GRID_PROP_MASK_GME,
  GRID_PROP_XVALS,
  GRID_PROP_YVALS,
  GRID_PROP_AREA,
  GRID_PROP_XBOUNDS,
  GRID_PROP_YBOUNDS,
};

struct gridVirtTable
{
  void (*destroy)(grid_t *gridptr);
  grid_t *(*copy)(grid_t *gridptr);
  void (*copyScalarFields)(grid_t *gridptrOrig, grid_t *gridptrDup);
  void (*copyArrayFields)(grid_t *gridptrOrig, grid_t *gridptrDup);
  void (*defXVals)(grid_t *gridptr, const double *xvals);
  void (*defYVals)(grid_t *gridptr, const double *yvals);
  void (*defMask)(grid_t *gridptr, const int *mask);
  void (*defMaskGME)(grid_t *gridptr, const int *mask);
  void (*defXBounds)(grid_t *gridptr, const double *xbounds);
  void (*defYBounds)(grid_t *gridptr, const double *ybounds);
  void (*defArea)(grid_t *gridptr, const double *area);
  double (*inqXVal)(grid_t *gridptr, SizeType index);
  double (*inqYVal)(grid_t *gridptr, SizeType index);
  SizeType (*inqXVals)(grid_t *gridptr, double *xvals);
  SizeType (*inqXValsPart)(grid_t *gridptr, int start, SizeType length, double *xvals);
  SizeType (*inqYVals)(grid_t *gridptr, double *yvals);
  SizeType (*inqYValsPart)(grid_t *gridptr, int start, SizeType length, double *yvals);
  const double *(*inqXValsPtr)(grid_t *gridptr);
  const double *(*inqYValsPtr)(grid_t *gridptr);
#ifndef USE_MPI
  int (*inqXIsc)(grid_t *gridptr);
  int (*inqYIsc)(grid_t *gridptr);
  SizeType (*inqXCvals)(grid_t *gridptr, char **xcvals);
  SizeType (*inqYCvals)(grid_t *gridptr, char **ycvals);
  const char **(*inqXCvalsPtr)(grid_t *gridptr);
  const char **(*inqYCvalsPtr)(grid_t *gridptr);
#endif
  double (*inqXInc)(grid_t *gridptr);
  double (*inqYInc)(grid_t *gridptr);
  // return true if for both grids, any xval and all yval differ
  bool (*compareXYFull)(grid_t *gridRef, grid_t *gridTest);
  // return if for both grids, x[0], y[0], x[size-1] and y[size-1] are respectively equal
  bool (*compareXYAO)(grid_t *gridRef, grid_t *gridTest);
  void (*inqArea)(grid_t *gridptr, double *area);
  const double *(*inqAreaPtr)(grid_t *gridptr);
  /* return 1 if inq property is set */
  int (*inqPropPresence)(grid_t *gridptr, enum gridPropInq inq);
  SizeType (*inqMask)(grid_t *gridptr, int *mask);
  int (*inqMaskGME)(grid_t *gridptr, int *mask_gme);
  SizeType (*inqXBounds)(grid_t *gridptr, double *xbounds);
  SizeType (*inqYBounds)(grid_t *gridptr, double *ybounds);
  const double *(*inqXBoundsPtr)(grid_t *gridptr);
  const double *(*inqYBoundsPtr)(grid_t *gridptr);
  int txCode;
  int (*getPackSize)(grid_t *gridptr, void *context);
  int (*getPackSizeScalars)(grid_t *gridptr, void *context);
  int (*getPackSizeArrays)(grid_t *gridptr, void *context);
  int (*unpack)(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id);
  grid_t *(*unpackScalars)(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context,
                           int force_id, int *memberMaskP);
  void (*unpackArrays)(grid_t *gridptr, int memberMask, char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos,
                       int originNamespace, void *context);
  void (*pack)(grid_t *gridptr, void *packBuffer, int packBufferSize, int *packBufferPos, void *context);
  /* return member mask */
  int (*packScalars)(grid_t *gridptr, void *packBuffer, int packBufferSize, int *packBufferPos, void *context);
  void (*packArrays)(grid_t *gridptr, int memberMask, void *packBuffer, int packBufferSize, int *packBufferPos, void *context);
};

struct gridaxis_t
{
  size_t size;  // number of values
  short flag;   // 0: undefined 1:vals 2:first+inc
  double first, last, inc;
  double *vals;
  double *bounds;
  cdi_keys_t keys;
#ifndef USE_MPI
  int clength;
  char **cvals;
#endif
};

// GME Grid
struct grid_gme_t
{
  int nd, ni, ni2, ni3;  // parameter for GRID_GME
};

struct grid_t
{
  char *name;
  int self;
  size_t size;
  int type;      // grid type
  int datatype;  // grid data type (used only internal in gridComplete())
  int proj;      // grid projection
  int projtype;  // grid projection type
  mask_t *mask;
  mask_t *mask_gme;
  double *area;
  struct grid_gme_t gme;
  int trunc;  // parameter for GRID_SPECTRAL
  int nvertex;
  int *reducedPoints;
  int reducedPointsSize;
  int np;                // number of parallels between a pole and the equator
  signed char isCyclic;  // three possible states:
                         // -1 if unknown,
                         //  0 if found not cyclic, or
                         //  1 for global cyclic grids
  bool lcomplex;
  bool hasdims;
  struct gridaxis_t x;
  struct gridaxis_t y;
  const struct gridVirtTable *vtable;
  cdi_keys_t keys;
  cdi_atts_t atts;
  void *extraData;
};

void grid_init(grid_t *gridptr);
void cdiGridTypeInit(grid_t *gridptr, int gridtype, size_t size);
void grid_free(grid_t *gridptr);
grid_t *grid_to_pointer(int gridID);
extern const struct gridVirtTable cdiGridVtable;

unsigned cdiGridCount(void);

void gridVerifyProj(int gridID);

double gridInqXincInMeter(int gridID);
double gridInqYincInMeter(int gridID);

// const double *gridInqXvalsPtr(int gridID);
// const double *gridInqYvalsPtr(int gridID);

const char **gridInqXCvalsPtr(int gridID);
const char **gridInqYCvalsPtr(int gridID);

// const double *gridInqXboundsPtr(int gridID);
// const double *gridInqYboundsPtr(int gridID);
const double *gridInqAreaPtr(int gridID);

int gridInqPropPresence(int gridID, enum gridPropInq inq);

int gridGenerate(const grid_t *grid);

// int gridIsEqual(int gridID1, int gridID2);

void cdiGridGetIndexList(unsigned, int *);

int gridUnpack(char *unpackBuffer, int unpackBufferSize, int *unpackBufferPos, int originNamespace, void *context, int force_id);

/* apply func to each grid */
enum cdiApplyRet cdiGridApply(enum cdiApplyRet (*func)(int id, void *res, void *data), void *data);

struct addIfNewRes
{
  int Id;
  int isNew;
};

struct addIfNewRes cdiVlistAddGridIfNew(int vlistID, grid_t *grid, int mode);

int gridVerifyProjParamsLCC(struct CDI_GridProjParams *gpp);
int gridVerifyProjParamsSTERE(struct CDI_GridProjParams *gpp);
int gridVerifyProjParamsHEALPIX(struct CDI_GridProjParams *gpp);

bool isGaussianLatitudes(size_t nlats, const double *latitudes);
void gaussianLatitudes(size_t nlats, double *latitudes, double *weights);

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
