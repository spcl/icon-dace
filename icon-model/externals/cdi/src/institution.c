#include <assert.h>
#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "namespace.h"
#include "serialize.h"
#include "institution.h"

typedef struct
{
  int self;
  int center;
  int subcenter;
  char *name;
  char *longname;
} institute_t;

static int instituteCompareKernel(institute_t *ip1, institute_t *ip2);
static void instituteDestroyP(institute_t *instituteptr);
static void institutePrintP(institute_t *instituteptr, FILE *fp);
static int instituteGetPackSize(institute_t *instituteptr, void *context);
static void institutePackP(void *instituteptr, void *buf, int size, int *position, void *context);
static int instituteTxCode(void *instituteptr);

static const resOps instituteOps = { (int (*)(void *, void *)) instituteCompareKernel,
                                     (void (*)(void *)) instituteDestroyP,
                                     (void (*)(void *, FILE *)) institutePrintP,
                                     (int (*)(void *, void *)) instituteGetPackSize,
                                     institutePackP,
                                     instituteTxCode };

static void
instituteDefaultValue(institute_t *instituteptr)
{
  instituteptr->self = CDI_UNDEFID;
  instituteptr->center = CDI_UNDEFID;
  instituteptr->subcenter = CDI_UNDEFID;
  instituteptr->name = NULL;
  instituteptr->longname = NULL;
}

void
instituteDefaultEntries(void)
{
  // clang-format off
  cdiResH resH[]
    = { institutDef( 98,   0, "ECMWF",     "European Centre for Medium-Range Weather Forecasts"),
        institutDef(252,   1, "MPIMET",    "Max Planck Institute for Meteorology"),
        institutDef( 98, 232, "MPIMET",    "Max Planck Institute for Meteorology"),
        institutDef( 98, 255, "MPIMET",    "Max-Planck-Institute for Meteorology"),
        institutDef( 78, 255, "DWD",       "Deutscher Wetterdienst"),
        institutDef( 78,   0, "DWD",       "Deutscher Wetterdienst"),
        institutDef(215, 255, "MCH",       "MeteoSwiss"),
        institutDef(  7,   0, "NCEP",      "National Centers for Environmental Prediction"),
        institutDef(  7,   1, "NCEP",      "National Centers for Environmental Prediction"),
        institutDef( 60,   0, "NCAR",      "National Center for Atmospheric Research"),
        institutDef( 74,   0, "METOFFICE", "U.K. Met Office"),
        institutDef( 97,   0, "ESA",       "European Space Agency"),
        institutDef( 99,   0, "KNMI",      "Royal Netherlands Meteorological Institute"),
        institutDef( 80,   0, "CNMC",      "Reparto per la Meteorologia, Rome (REMET)"),
        // institutDef(  0,   0, "IPSL", "IPSL (Institut Pierre Simon Laplace, Paris, France)");
  };
  // clang-format on

  const size_t n = sizeof(resH) / sizeof(*resH);
  for (size_t i = 0; i < n; i++) reshSetStatus(resH[i], &instituteOps, RESH_IN_USE);
}

static int
instituteCompareKernel(institute_t *ip1, institute_t *ip2)
{
  int differ = 0;

  if (ip1->name)
    {
      if (ip1->center > 0 && ip2->center != ip1->center) differ = 1;
      if (ip1->subcenter > 0 && ip2->subcenter != ip1->subcenter) differ = 1;

      if (!differ)
        {
          if (ip2->name)
            {
              const size_t len1 = strlen(ip1->name);
              const size_t len2 = strlen(ip2->name);
              if ((len1 != len2) || memcmp(ip2->name, ip1->name, len2)) differ = 1;
            }
        }
    }
  else if (ip1->longname)
    {
      if (ip2->longname)
        {
          const size_t len1 = strlen(ip1->longname);
          const size_t len2 = strlen(ip2->longname);
          if ((len1 != len2) || memcmp(ip2->longname, ip1->longname, len2)) differ = 1;
        }
    }
  else
    {
      if (!(ip2->center == ip1->center && ip2->subcenter == ip1->subcenter)) differ = 1;
      if (ip1->subcenter > 0 && ip1->subcenter != 255 && ip2->subcenter != ip1->subcenter) differ = 1;
    }

  return differ;
}

struct instLoc
{
  institute_t *ip;
  int id;
};

static enum cdiApplyRet
findInstitute(int id, void *res, void *data)
{
  institute_t *ip1 = ((struct instLoc *) data)->ip;
  institute_t *ip2 = (institute_t *) res;
  if (!instituteCompareKernel(ip1, ip2))
    {
      ((struct instLoc *) data)->id = id;
      return CDI_APPLY_STOP;
    }
  else
    return CDI_APPLY_GO_ON;
}

int
institutInq(int center, int subcenter, const char *name, const char *longname)
{
  institute_t ip_ref;
  ip_ref.self = CDI_UNDEFID;
  ip_ref.center = center;
  ip_ref.subcenter = subcenter;
  ip_ref.name = (name && name[0]) ? (char *) name : NULL;
  ip_ref.longname = (longname && longname[0]) ? (char *) longname : NULL;

  struct instLoc state = { .ip = &ip_ref, .id = CDI_UNDEFID };
  cdiResHFilterApply(&instituteOps, findInstitute, &state);

  return state.id;
}

static institute_t *
instituteNewEntry(cdiResH resH, int center, int subcenter, const char *name, const char *longname)
{
  institute_t *instituteptr = (institute_t *) Malloc(sizeof(institute_t));
  instituteDefaultValue(instituteptr);
  if (resH == CDI_UNDEFID)
    instituteptr->self = reshPut(instituteptr, &instituteOps);
  else
    {
      instituteptr->self = resH;
      reshReplace(resH, instituteptr, &instituteOps);
    }
  instituteptr->center = center;
  instituteptr->subcenter = subcenter;
  if (name && *name) instituteptr->name = strdup(name);
  if (longname && *longname) instituteptr->longname = strdup(longname);
  return instituteptr;
}

int
institutDef(int center, int subcenter, const char *name, const char *longname)
{
  institute_t *instituteptr = instituteNewEntry(CDI_UNDEFID, center, subcenter, name, longname);
  return instituteptr->self;
}

int
institutInqCenter(int instID)
{
  return instID != CDI_UNDEFID ? ((institute_t *) (reshGetVal(instID, &instituteOps)))->center : CDI_UNDEFID;
}

int
institutInqSubcenter(int instID)
{
  return instID != CDI_UNDEFID ? ((institute_t *) (reshGetVal(instID, &instituteOps)))->subcenter : CDI_UNDEFID;
}

const char *
institutInqNamePtr(int instID)
{
  return instID != CDI_UNDEFID ? ((institute_t *) (reshGetVal(instID, &instituteOps)))->name : NULL;
}

const char *
institutInqLongnamePtr(int instID)
{
  return instID != CDI_UNDEFID ? ((institute_t *) (reshGetVal(instID, &instituteOps)))->longname : NULL;
}

int
institutInqNumber(void)
{
  int instNum = (int) (reshCountType(&instituteOps));
  return instNum;
}

static void
instituteDestroyP(institute_t *instituteptr)
{
  xassert(instituteptr);
  Free(instituteptr->name);
  Free(instituteptr->longname);
  Free(instituteptr);
}

static void
institutePrintP(institute_t *ip, FILE *fp)
{
  if (ip)
    fprintf(fp,
            "#\n"
            "# instituteID %d\n"
            "#\n"
            "self          = %d\n"
            "center        = %d\n"
            "subcenter     = %d\n"
            "name          = %s\n"
            "longname      = %s\n",
            ip->self, ip->self, ip->center, ip->subcenter, ip->name ? ip->name : "NN", ip->longname ? ip->longname : "NN");
}

static int
instituteTxCode(void *instituteptr)
{
  (void) instituteptr;
  return INSTITUTE;
}

enum
{
  INSTITUTE_PACK_INT_SELF,
  INSTITUTE_PACK_INT_CENTER,
  INSTITUTE_PACK_INT_SUBCENTER,
  INSTITUTE_PACK_INT_NAMELEN,
  INSTITUTE_PACK_INT_LNAMELEN,
  institute_nints,
};

static int
instituteGetPackSize(institute_t *ip, void *context)
{
  size_t namelen = strlen(ip->name), longnamelen = strlen(ip->longname);
  xassert(namelen < INT_MAX && longnamelen < INT_MAX);
  size_t txsize = (size_t) serializeGetSize(institute_nints, CDI_DATATYPE_INT, context)
                  + (size_t) serializeGetSize((int) namelen + 1, CDI_DATATYPE_TXT, context)
                  + (size_t) serializeGetSize((int) longnamelen + 1, CDI_DATATYPE_TXT, context);
  xassert(txsize <= INT_MAX);
  return (int) txsize;
}

static void
institutePackP(void *instituteptr, void *buf, int size, int *position, void *context)
{
  institute_t *p = (institute_t *) instituteptr;
  int tempbuf[institute_nints];
  tempbuf[INSTITUTE_PACK_INT_SELF] = p->self;
  tempbuf[INSTITUTE_PACK_INT_CENTER] = p->center;
  tempbuf[INSTITUTE_PACK_INT_SUBCENTER] = p->subcenter;
  tempbuf[INSTITUTE_PACK_INT_NAMELEN] = (int) strlen(p->name) + 1;
  tempbuf[INSTITUTE_PACK_INT_LNAMELEN] = (int) strlen(p->longname) + 1;
  serializePack(tempbuf, institute_nints, CDI_DATATYPE_INT, buf, size, position, context);
  serializePack(p->name, tempbuf[INSTITUTE_PACK_INT_NAMELEN], CDI_DATATYPE_TXT, buf, size, position, context);
  serializePack(p->longname, tempbuf[INSTITUTE_PACK_INT_LNAMELEN], CDI_DATATYPE_TXT, buf, size, position, context);
}

int
instituteUnpack(void *buf, int size, int *position, int originNamespace, void *context, int force_id)
{
#define adaptKey(key) (namespaceAdaptKey((key), originNamespace))
  int tempbuf[institute_nints];
  int instituteID;
  serializeUnpack(buf, size, position, tempbuf, institute_nints, CDI_DATATYPE_INT, context);
  char *name = (char *) Malloc((size_t) tempbuf[INSTITUTE_PACK_INT_NAMELEN] + (size_t) tempbuf[INSTITUTE_PACK_INT_LNAMELEN]),
       *longname = name + tempbuf[INSTITUTE_PACK_INT_NAMELEN];
  serializeUnpack(buf, size, position, name, tempbuf[INSTITUTE_PACK_INT_NAMELEN], CDI_DATATYPE_TXT, context);
  serializeUnpack(buf, size, position, longname, tempbuf[INSTITUTE_PACK_INT_LNAMELEN], CDI_DATATYPE_TXT, context);
  int targetID = force_id ? adaptKey(tempbuf[INSTITUTE_PACK_INT_SELF]) : CDI_UNDEFID;
  institute_t *ip
      = instituteNewEntry(targetID, tempbuf[INSTITUTE_PACK_INT_CENTER], tempbuf[INSTITUTE_PACK_INT_SUBCENTER], name, longname);
  instituteID = ip->self;
  xassert(!force_id || instituteID == targetID);
  Free(name);
  reshSetStatus(instituteID, &instituteOps, reshGetStatus(instituteID, &instituteOps) & ~RESH_SYNC_BIT);
#undef adaptKey
  return instituteID;
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
