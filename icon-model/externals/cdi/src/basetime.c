#include <stddef.h>  // for NULL
#include "cdi.h"
#include "basetime.h"
#include "error.h"

void
basetimeInit(basetime_t *basetime)
{
  if (basetime == NULL) Error("Internal problem! Basetime not allocated.");

  if (basetime)
    {
      basetime->ncvarid = CDI_UNDEFID;
      basetime->ncdimid = CDI_UNDEFID;
      basetime->ncvarboundsid = CDI_UNDEFID;
      basetime->leadtimeid = CDI_UNDEFID;
      basetime->hasUnits = false;
      basetime->isWRF = false;
    }
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
