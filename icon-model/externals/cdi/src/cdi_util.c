#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <unistd.h>

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "cdi.h"
#include "cdi_int.h"

void
cdiDecodeParam(int param, int *pnum, int *pcat, int *pdis)
{
  unsigned uparam = (unsigned) param;

  *pdis = (int) (0xffU & uparam);
  *pcat = (int) (0xffU & (uparam >> 8));
  unsigned upnum = 0xffffU & (uparam >> 16);
  if (upnum > 0x7fffU) upnum = 0x8000U - upnum;
  *pnum = (int) upnum;
}

int
cdiEncodeParam(int pnum, int pcat, int pdis)
{
  if (pcat < 0 || pcat > 255) pcat = 255;
  if (pdis < 0 || pdis > 255) pdis = 255;

  unsigned upnum = (unsigned) pnum;
  if (pnum < 0) upnum = (unsigned) (0x8000 - pnum);

  unsigned uparam = (upnum << 16) | (((unsigned) pcat) << 8) | (unsigned) pdis;

  return (int) uparam;
}

void
cdiParamToString(int param, char *paramstr, int maxlen)
{
  int dis, cat, num;
  cdiDecodeParam(param, &num, &cat, &dis);

  size_t umaxlen = maxlen >= 0 ? (unsigned) maxlen : 0U;
  int len;
  if (dis == 255 && (cat == 255 || cat == 0))
    len = snprintf(paramstr, umaxlen, "%d", num);
  else if (dis == 255)
    len = snprintf(paramstr, umaxlen, "%d.%d", num, cat);
  else
    len = snprintf(paramstr, umaxlen, "%d.%d.%d", num, cat, dis);

  if (len >= maxlen || len < 0) fprintf(stderr, "Internal problem (%s): size of input string is too small!\n", __func__);
}

const char *
cdiUnitNamePtr(int cdi_unit)
{
  const char *cdiUnits[] = {
    /*  0 */ "undefined",
    /*  1 */ "Pa",
    /*  2 */ "hPa",
    /*  3 */ "mm",
    /*  4 */ "cm",
    /*  5 */ "dm",
    /*  6 */ "m",
  };
  enum
  {
    numUnits = (int) (sizeof(cdiUnits) / sizeof(char *))
  };
  const char *name = (cdi_unit > 0 && cdi_unit < numUnits) ? cdiUnits[cdi_unit] : NULL;

  return name;
}

size_t
cdiGetPageSize(bool largePageAlign)
{
  long pagesize = -1L;
#if HAVE_DECL__SC_LARGE_PAGESIZE || HAVE_DECL__SC_PAGE_SIZE || HAVE_DECL__SC_PAGESIZE
  bool nameAssigned = false;
  int name = 0;
#if HAVE_DECL__SC_LARGE_PAGESIZE
  if (largePageAlign)
    {
      name = _SC_LARGE_PAGESIZE;
      nameAssigned = true;
    }
  else
#else
  (void) largePageAlign;
#endif
    {
#if HAVE_DECL__SC_PAGESIZE || HAVE_DECL__SC_PAGE_SIZE
      name =
#if HAVE_DECL__SC_PAGESIZE
          _SC_PAGESIZE;
#else
          _SC_PAGE_SIZE;
#endif
      nameAssigned = true;
#endif
    }
  if (nameAssigned) pagesize = sysconf(name);
#endif
  if (pagesize == -1L)
    pagesize =
#if HAVE_DECL_PAGESIZE
        PAGESIZE;
#elif HAVE_DECL_PAGE_SIZE
        PAGE_SIZE;
#else
        commonPageSize;
#endif
  return (size_t) pagesize;
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
