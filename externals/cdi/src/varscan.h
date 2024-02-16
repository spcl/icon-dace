#ifndef VARSCAN_H
#define VARSCAN_H

#ifndef GRID_H
#include "grid.h"
#endif

void varAddRecord(int recID, int param, int gridID, int zaxistype, int lbounds, int level1, int level2, int level_sf,
                  int level_unit, int prec, int *pvarID, int *plevelID, int tsteptype, int ltype1, int ltype2, const char *name,
                  const VarScanKeys *scanKeys, const var_tile_t *tiles, int *tile_index);

void varDefVCT(size_t vctsize, double *vctptr);
void varDefZAxisReference(int nlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE]);

int varDefZaxis(int vlistID, int zaxistype, int nlevels, const double *levels, const char **cvals, size_t clength, bool lbounds,
                const double *levels1, const double *levels2, int vctsize, const double *vct, char *name, const char *longname,
                const char *units, int prec, int mode, int ltype1, int ltype2);

void varDefMissval(int varID, double missval);
void varDefCompType(int varID, int comptype);
void varDefCompLevel(int varID, int complevel);
void varDefInst(int varID, int instID);
int varInqInst(int varID);
void varDefModel(int varID, int modelID);
int varInqModel(int varID);
void varDefTable(int varID, int tableID);
int varInqTable(int varID);

void varDefKeyInt(int varID, int key, int value);
void varDefKeyBytes(int varID, int key, const unsigned char *bytes, int length);
void varDefKeyString(int varID, int key, const char *string);

void varDefOptGribInt(int varID, int tile_index, long lval, const char *keyword);
void varDefOptGribDbl(int varID, int tile_index, double dval, const char *keyword);
int varOptGribNentries(int varID);

bool zaxis_compare(int zaxisID, int zaxistype, int nlevels, const double *levels, const double *lbounds, const double *ubounds,
                   const char *longname, const char *units, int ltype1, int ltype2);

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
