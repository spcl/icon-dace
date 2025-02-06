#ifndef TABLEPAR_H
#define TABLEPAR_H

enum
{
  TABLE_DUP_NAME = 1 << 0,
  TABLE_DUP_LONGNAME = 1 << 1,
  TABLE_DUP_UNITS = 1 << 2,
};

typedef struct
{
  int id;                // Parameter number (GRIB)
  int ltype;             // Level type (GRIB)
  int dupflags;          // keep track of which attributes got strdup'ed
  const char *name;      // Parameter name
  const char *longname;  // Parameter long name
  const char *units;     // Parameter units
} param_type;

void tableLink(int tableID, const param_type *pars, int npars);
int tableDef(int modelID, int tablegribID, const char *tablename);

int tableInqParCode(int tableID, char *name, int *code);

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
