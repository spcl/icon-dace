#include <stddef.h>
#include <string.h>
#include <ctype.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"

#include "tablepar.h"
#include "table.h"

#define MAX_TABLE 256
#define MAX_PARS 1024

typedef struct
{
  bool used;
  int npars;
  int modelID;
  int number;
  char *name;
  param_type *pars;
} paramtab_type;

static paramtab_type parTable[MAX_TABLE];
static int parTableSize = MAX_TABLE;
static int parTableNum = 0;
static int ParTableInit = 0;

static char *tablePath = NULL;

static void tableDefModelID(int tableID, int modelID);
static void tableDefNum(int tableID, int tablenum);

static void
tableDefEntry(int tableID, int id, int ltype, const char *name, const char *longname, const char *units)
{
  if (tableID >= 0 && tableID < MAX_TABLE && parTable[tableID].used)
    {
    }
  else
    Error("Invalid table ID %d", tableID);

  int item = parTable[tableID].npars++;
  parTable[tableID].pars[item].id = id;
  parTable[tableID].pars[item].ltype = ltype;
  parTable[tableID].pars[item].dupflags = 0;
  parTable[tableID].pars[item].name = NULL;
  parTable[tableID].pars[item].longname = NULL;
  parTable[tableID].pars[item].units = NULL;

  if (name && name[0])
    {
      parTable[tableID].pars[item].name = strdupx(name);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_NAME;
    }
  if (longname && longname[0])
    {
      parTable[tableID].pars[item].longname = strdupx(longname);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_LONGNAME;
    }
  if (units && units[0])
    {
      parTable[tableID].pars[item].units = strdupx(units);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_UNITS;
    }
}

void
tableLink(int tableID, const param_type *pars, int npars)
{
  for (int item = 0; item < npars; item++)
    {
      parTable[tableID].pars[item].id = pars[item].id;
      parTable[tableID].pars[item].ltype = pars[item].ltype;
      parTable[tableID].pars[item].dupflags = 0;
      parTable[tableID].pars[item].name = pars[item].name;
      parTable[tableID].pars[item].longname = pars[item].longname;
      parTable[tableID].pars[item].units = pars[item].units;
    }

  parTable[tableID].npars = npars;
}

static void
parTableInitEntry(int tableID)
{
  parTable[tableID].used = false;
  parTable[tableID].pars = NULL;
  parTable[tableID].npars = 0;
  parTable[tableID].modelID = CDI_UNDEFID;
  parTable[tableID].number = CDI_UNDEFID;
  parTable[tableID].name = NULL;
}

static void
tableGetPath(void)
{
  char *path = getenv("TABLEPATH");
  if (path) tablePath = strdupx(path);
  // printf("tablePath = %s\n", tablePath);
}

static void
parTableFinalize(void)
{
  for (int tableID = 0; tableID < MAX_TABLE; ++tableID)
    if (parTable[tableID].used)
      {
        int npars = parTable[tableID].npars;
        for (int item = 0; item < npars; ++item)
          {
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_NAME) Free((void *) parTable[tableID].pars[item].name);
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_LONGNAME) Free((void *) parTable[tableID].pars[item].longname);
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_UNITS) Free((void *) parTable[tableID].pars[item].units);
          }
        Free(parTable[tableID].pars);
        Free(parTable[tableID].name);
      }
}

static void
parTableInit(void)
{
  ParTableInit = 1;

  atexit(parTableFinalize);
  if (cdiPartabIntern) tableDefault();

  tableGetPath();
}

static int
tableNewEntry(void)
{
  int tableID = 0;

  static int init = 0;
  if (!init)
    {
      for (tableID = 0; tableID < parTableSize; tableID++) parTableInitEntry(tableID);
      init = 1;
    }

  // Look for a free slot in parTable.
  for (tableID = 0; tableID < parTableSize; tableID++)
    {
      if (!parTable[tableID].used) break;
    }

  if (tableID == parTableSize) Error("no more entries!");

  parTable[tableID].used = true;
  parTableNum++;

  return tableID;
}

static int
decodeForm1(char *pline, char *name, char *longname, char *units)
{
  char *pstart, *pend;

  // FIXME: parse success isn't verified
  /* long level =  */ strtol(pline, &pline, 10);
  while (isspace((int) *pline)) pline++;

  pstart = pline;
  while (!(isspace((int) *pline) || *pline == 0)) pline++;
  size_t len = (size_t) (pline - pstart);
  if (len > 0)
    {
      memcpy(name, pstart, len);
      name[len] = 0;
    }
  else
    return 0;

  if (pline[0] == 0) return 0;

  // Format 1 : code name add mult longname [units]
  // FIXME: successful parse isn't verified
  /* double add  =  */ strtod(pline, &pline);
  // FIXME: successful parse isn't verified
  /* double mult =  */ strtod(pline, &pline);

  while (isspace((int) *pline)) pline++;

  len = strlen(pline);
  if (len > 0)
    {
      pstart = pline;
      pend = strrchr(pline, '[');
      if (pend == pstart)
        len = 0;
      else
        {
          if (pend)
            pend--;
          else
            pend = pstart + len;
          while (isspace((int) *pend)) pend--;
          len = (size_t) (pend - pstart + 1);
        }
      if (len > 0)
        {
          memcpy(longname, pstart, len);
          longname[len] = 0;
        }
      pstart = strrchr(pline, '[');
      if (pstart)
        {
          pstart++;
          while (isspace((int) *pstart)) pstart++;
          pend = strchr(pstart, ']');
          if (!pend) return 0;
          pend--;
          while (isspace((int) *pend)) pend--;
          len = (size_t) (pend - pstart + 1);
          if (len > 0)
            {
              memcpy(units, pstart, len);
              units[len] = 0;
            }
        }
    }

  return 0;
}

static int
decodeForm2(char *pline, char *name, char *longname, char *units)
{
  // Format 2 : code | name | longname | units
  char *pend;

  pline = strchr(pline, '|');
  pline++;

  while (isspace((int) *pline)) pline++;
  if (*pline != '|')
    {
      pend = strchr(pline, '|');
      if (!pend)
        {
          pend = pline;
          while (!isspace((int) *pend)) pend++;
          size_t len = (size_t) (pend - pline);
          if (len > 0)
            {
              memcpy(name, pline, len);
              name[len] = 0;
            }
          return 0;
        }
      else
        {
          pend--;
          while (isspace((int) *pend)) pend--;
          size_t len = (size_t) (pend - pline + 1);
          if (len > 0)
            {
              memcpy(name, pline, len);
              name[len] = 0;
            }
        }
    }
  else
    name[0] = '\0';

  pline = strchr(pline, '|');
  pline++;
  while (isspace((int) *pline)) pline++;
  pend = strchr(pline, '|');
  if (!pend) pend = strchr(pline, 0);
  pend--;
  while (isspace((int) *pend)) pend--;
  {
    size_t len = (size_t) (pend - pline + 1);
    if (len > 0)
      {
        memcpy(longname, pline, len);
        longname[len] = 0;
      }
  }

  pline = strchr(pline, '|');
  if (pline)
    {
      pline++;
      while (isspace((int) *pline)) pline++;
      pend = strchr(pline, '|');
      if (!pend) pend = strchr(pline, 0);
      pend--;
      while (isspace((int) *pend)) pend--;
      ptrdiff_t len = pend - pline + 1;
      if (len < 0) len = 0;
      memcpy(units, pline, (size_t) len);
      units[len] = 0;
    }

  return 0;
}

int
tableRead(const char *tablefile)
{
  char line[1024], *pline;
  char name[256], longname[256], units[256];
  int tableID = CDI_UNDEFID;

  FILE *tablefp = fopen(tablefile, "r");
  if (tablefp == NULL) return tableID;

  const char *tablename = strrchr(tablefile, '/');
  if (tablename == 0)
    tablename = tablefile;
  else
    tablename++;

  tableID = tableDef(-1, 0, tablename);

  while (fgets(line, 1023, tablefp))
    {
      size_t len = strlen(line);
      if (line[len - 1] == '\n') line[len - 1] = '\0';
      name[0] = 0;
      longname[0] = 0;
      units[0] = 0;
      if (line[0] == '#') continue;
      pline = line;

      len = strlen(pline);
      if (len < 4) continue;
      while (isspace((int) *pline)) pline++;
      int id = atoi(pline);
      // if ( id > 255 ) id -= 256;
      if (id == 0) continue;

      while (isdigit((int) *pline)) pline++;

      int ltype = CDI_UNDEFID;
      if (*pline == ';' || *pline == ':')
        {
          pline++;
          ltype = atoi(pline);
          while (isdigit((int) *pline)) pline++;

          if (*pline == ';' || *pline == ':')
            {
              pline++;
              while (isdigit((int) *pline)) pline++;
            }
        }

      while (isdigit((int) *pline)) pline++;

      int err = (strchr(pline, '|')) ? decodeForm2(pline, name, longname, units) : decodeForm1(pline, name, longname, units);
      if (err) continue;

      if (name[0] == 0) snprintf(name, sizeof(name), "var%d", id);

      tableDefEntry(tableID, id, ltype, name, longname, units);
    }

  return tableID;
}

static int
tableFromEnv(int modelID, int tablenum)
{
  char tablename[256] = { '\0' };
  size_t tablenameLen = 0;
  int instID;

  const char *name2Use;
  {
    const char *modelName, *instName;
    if ((modelName = modelInqNamePtr(modelID)))
      name2Use = modelName;
    else if ((instID = modelInqInstitut(modelID)) != CDI_UNDEFID && (instName = institutInqNamePtr(instID)))
      name2Use = instName;
    else
      return CDI_UNDEFID;
  }
  tablenameLen = strlen(name2Use);
  memcpy(tablename, name2Use, tablenameLen);
  if (tablenum) tablenameLen += (size_t) (sprintf(tablename + tablenameLen, "_%03d", tablenum));
  size_t lenp = 0, lenf = tablenameLen;
  if (tablePath) lenp = strlen(tablePath);
  // if (tablePath) printf("tablePath = %s\n", tablePath);
  // if (tablename) printf("tableName = %s\n", tablename);
  char *tablefile = (char *) Malloc(lenp + lenf + 3);
  if (tablePath)
    {
      strcpy(tablefile, tablePath);
      strcat(tablefile, "/");
    }
  else
    tablefile[0] = '\0';
  strcat(tablefile, tablename);
  // if (tablefile) printf("tableFile = %s\n", tablefile);

  int tableID = tableRead(tablefile);
  if (tableID != CDI_UNDEFID)
    {
      tableDefModelID(tableID, modelID);
      tableDefNum(tableID, tablenum);
    }
  // printf("tableID = %d %s\n", tableID, tablefile);
  Free(tablefile);

  return tableID;
}

int
tableInq(int modelID, int tablenum, const char *tablename)
{
  int tableID = CDI_UNDEFID;
  int modelID2 = CDI_UNDEFID;
  char tablefile[256] = { '\0' };

  if (!ParTableInit) parTableInit();

  if (tablename)
    {
      strcpy(tablefile, tablename);
      /*
      printf("tableInq: tablefile = >%s<\n", tablefile);
      */
      /* search for internal table */
      for (tableID = 0; tableID < MAX_TABLE; tableID++)
        {
          if (parTable[tableID].used && parTable[tableID].name)
            {
              /* len = strlen(parTable[tableID].name); */
              size_t len = strlen(tablename);
              if (memcmp(parTable[tableID].name, tablename, len) == 0) break;
            }
        }
      if (tableID == MAX_TABLE) tableID = CDI_UNDEFID;
      if (CDI_Debug) Message("tableID = %d tablename = %s", tableID, tablename);
    }
  else
    {
      for (tableID = 0; tableID < MAX_TABLE; tableID++)
        {
          if (parTable[tableID].used)
            {
              if (parTable[tableID].modelID == modelID && parTable[tableID].number == tablenum) break;
            }
        }

      if (tableID == MAX_TABLE) tableID = CDI_UNDEFID;

      if (tableID == CDI_UNDEFID)
        {
          if (modelID != CDI_UNDEFID)
            {
              const char *modelName;
              if ((modelName = modelInqNamePtr(modelID)))
                {
                  strcpy(tablefile, modelName);
                  size_t len = strlen(tablefile);
                  for (size_t i = 0; i < len; i++)
                    if (tablefile[i] == '.') tablefile[i] = '\0';
                  modelID2 = modelInq(-1, 0, tablefile);
                }
            }
          if (modelID2 != CDI_UNDEFID)
            for (tableID = 0; tableID < MAX_TABLE; tableID++)
              {
                if (parTable[tableID].used)
                  {
                    if (parTable[tableID].modelID == modelID2 && parTable[tableID].number == tablenum) break;
                  }
              }
        }

      if (tableID == MAX_TABLE) tableID = CDI_UNDEFID;

      if (tableID == CDI_UNDEFID && modelID != CDI_UNDEFID) tableID = tableFromEnv(modelID, tablenum);

      if (CDI_Debug && tablename) Message("tableID = %d tablename = %s", tableID, tablename);
    }

  return tableID;
}

int
tableDef(int modelID, int tablenum, const char *tablename)
{
  int tableID = CDI_UNDEFID;

  if (!ParTableInit) parTableInit();
  // if (!(modelID == CDI_UNDEFID && tablenum == 0)) tableID = tableInq(modelID, tablenum, tablename);
  if (tableID == CDI_UNDEFID)
    {
      tableID = tableNewEntry();

      parTable[tableID].modelID = modelID;
      parTable[tableID].number = tablenum;
      if (tablename) parTable[tableID].name = strdupx(tablename);

      parTable[tableID].pars = (param_type *) Malloc(MAX_PARS * sizeof(param_type));
    }

  return tableID;
}

static void
tableDefModelID(int tableID, int modelID)
{
  parTable[tableID].modelID = modelID;
}

static void
tableDefNum(int tableID, int tablenum)
{
  parTable[tableID].number = tablenum;
}

int
tableInqNum(int tableID)
{
  int number = 0;

  if (tableID >= 0 && tableID < MAX_TABLE) number = parTable[tableID].number;

  return number;
}

int
tableInqModel(int tableID)
{
  int modelID = -1;

  if (tableID >= 0 && tableID < MAX_TABLE) modelID = parTable[tableID].modelID;

  return modelID;
}

static void
partabCheckID(int item)
{
  if (item < 0 || item >= parTableSize) Error("item %d undefined!", item);

  if (!parTable[item].name) Error("item %d name undefined!", item);
}

const char *
tableInqNamePtr(int tableID)
{
  const char *tablename = NULL;

  if (CDI_Debug) Message("tableID = %d", tableID);

  if (!ParTableInit) parTableInit();

  if (tableID >= 0 && tableID < parTableSize)
    if (parTable[tableID].name) tablename = parTable[tableID].name;

  return tablename;
}

static size_t
max_length(size_t maxlen, const char *cstring)
{
  if (cstring)
    {
      size_t len = strlen(cstring);
      if (len > maxlen) maxlen = len;
    }
  return maxlen;
}

void
tableWrite(const char *ptfile, int tableID)
{
  size_t maxname = 4, maxlname = 10, maxunits = 2;
  int instID = CDI_UNDEFID;
  int center = 0, subcenter = 0;
  const char *instnameptr = NULL, *modelnameptr = NULL;

  if (CDI_Debug) Message("write parameter table %d to %s", tableID, ptfile);

  if (tableID == CDI_UNDEFID)
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  FILE *ptfp = fopen(ptfile, "w");

  int npars = parTable[tableID].npars;
  for (int item = 0; item < npars; item++)
    {
      maxname = max_length(maxname, parTable[tableID].pars[item].name);
      maxlname = max_length(maxlname, parTable[tableID].pars[item].longname);
      maxunits = max_length(maxunits, parTable[tableID].pars[item].units);
    }

  int tablenum = tableInqNum(tableID);
  int modelID = parTable[tableID].modelID;
  if (modelID != CDI_UNDEFID)
    {
      modelnameptr = modelInqNamePtr(modelID);
      instID = modelInqInstitut(modelID);
    }
  if (instID != CDI_UNDEFID)
    {
      center = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);
      instnameptr = institutInqNamePtr(instID);
    }

  fprintf(ptfp, "# Parameter table\n");
  fprintf(ptfp, "#\n");
  if (tablenum) fprintf(ptfp, "# TABLE_ID=%d\n", tablenum);
  fprintf(ptfp, "# TABLE_NAME=%s\n", parTable[tableID].name);
  if (modelnameptr) fprintf(ptfp, "# TABLE_MODEL=%s\n", modelnameptr);
  if (instnameptr) fprintf(ptfp, "# TABLE_INSTITUT=%s\n", instnameptr);
  if (center) fprintf(ptfp, "# TABLE_CENTER=%d\n", center);
  if (subcenter) fprintf(ptfp, "# TABLE_SUBCENTER=%d\n", subcenter);
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id       = parameter ID\n");
  fprintf(ptfp, "# name     = variable name\n");
  fprintf(ptfp, "# title    = long name (description)\n");
  fprintf(ptfp, "# units    = variable units\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# The format of each record is:\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id | %-*s | %-*s | %-*s\n", (int) maxname, "name", (int) maxlname, "title", (int) maxunits, "units");

  for (int item = 0; item < npars; item++)
    {
      const char *name = parTable[tableID].pars[item].name, *longname = parTable[tableID].pars[item].longname,
                 *units = parTable[tableID].pars[item].units;
      if (name == NULL) name = " ";
      if (longname == NULL) longname = " ";
      if (units == NULL) units = " ";
      fprintf(ptfp, "%4d | %-*s | %-*s | %-*s\n", parTable[tableID].pars[item].id, (int) maxname, name, (int) maxlname, longname,
              (int) maxunits, units);
    }

  fclose(ptfp);
}

void
tableFWriteC(FILE *ptfp, int tableID)
{
  const char chelp[] = "";
  size_t maxname = 0, maxlname = 0, maxunits = 0;
  char tablename[256];

  if (tableID == CDI_UNDEFID)
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  int npars = parTable[tableID].npars;
  for (int item = 0; item < npars; item++)
    {
      maxname = max_length(maxname, parTable[tableID].pars[item].name);
      maxlname = max_length(maxlname, parTable[tableID].pars[item].longname);
      maxunits = max_length(maxunits, parTable[tableID].pars[item].units);
    }

  strncpy(tablename, parTable[tableID].name, sizeof(tablename) - 1);
  tablename[sizeof(tablename) - 1] = '\0';
  {
    size_t len = strlen(tablename);
    for (size_t i = 0; i < len; i++)
      if (tablename[i] == '.') tablename[i] = '_';
  }
  fprintf(ptfp, "static const param_type %s[] = {\n", tablename);

  for (int item = 0; item < npars; item++)
    {
      size_t len = strlen(parTable[tableID].pars[item].name),
             llen = parTable[tableID].pars[item].longname ? strlen(parTable[tableID].pars[item].longname) : 0,
             ulen = parTable[tableID].pars[item].units ? strlen(parTable[tableID].pars[item].units) : 0;
      fprintf(ptfp, "  {%4d, -1, 0, \"%s\", %-*s%c%s%s, %-*s%c%s%s %-*s},\n", parTable[tableID].pars[item].id,
              parTable[tableID].pars[item].name, (int) (maxname - len), chelp, llen ? '"' : ' ',
              llen ? parTable[tableID].pars[item].longname : "NULL", llen ? "\"" : "", (int) (maxlname - (llen ? llen : 3)), chelp,
              ulen ? '"' : ' ', ulen ? parTable[tableID].pars[item].units : "NULL", ulen ? "\"" : "",
              (int) (maxunits - (ulen ? ulen : 3)), chelp);
    }

  fprintf(ptfp, "};\n\n");
}

void
tableInqEntry(int tableID, int id, int ltype, char *name, char *longname, char *units)
{
  if (((tableID >= 0) & (tableID < MAX_TABLE)) | (tableID == CDI_UNDEFID))
    {
    }
  else
    Error("Invalid table ID %d", tableID);

  if (tableID != CDI_UNDEFID)
    {
      int npars = parTable[tableID].npars;
      for (int item = 0; item < npars; item++)
        {
          if (parTable[tableID].pars[item].id == id
              && (parTable[tableID].pars[item].ltype == -1 || ltype == -1 || parTable[tableID].pars[item].ltype == ltype))
            {
              if (name && parTable[tableID].pars[item].name) strcpy(name, parTable[tableID].pars[item].name);
              if (longname && parTable[tableID].pars[item].longname) strcpy(longname, parTable[tableID].pars[item].longname);
              if (units && parTable[tableID].pars[item].units) strcpy(units, parTable[tableID].pars[item].units);

              break;
            }
        }
    }
}

int
tableInqParCode(int tableID, char *varname, int *code)
{
  int err = 1;

  if (tableID != CDI_UNDEFID && varname != NULL)
    {
      int npars = parTable[tableID].npars;
      for (int item = 0; item < npars; item++)
        {
          if (parTable[tableID].pars[item].name && str_is_equal(parTable[tableID].pars[item].name, varname))
            {
              *code = parTable[tableID].pars[item].id;
              err = 0;
              break;
            }
        }
    }

  return err;
}

int
tableInqNumber(void)
{
  if (!ParTableInit) parTableInit();

  return parTableNum;
}
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
