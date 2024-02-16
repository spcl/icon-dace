#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h>
#endif

#include <stdio.h>

#include "cdi.h"
#include "cdi_int.h"
#include "gribapi.h"
#include "dmemory.h"

static char gribapi_libvers[64] = "";
#ifdef HAVE_LIBGRIB_API
static bool gribapi_libvers_init;
#endif

void
gribapiLibraryVersion(int *major_version, int *minor_version, int *revision_version)
{
#ifdef HAVE_LIBGRIB_API
  long version = grib_get_api_version();
  (*major_version) = (int) (version / 10000);
  (*minor_version) = (int) ((version - (*major_version) * 10000) / 100);
  (*revision_version) = (int) (version - (*major_version) * 10000 - (*minor_version) * 100);
#else
  (*major_version) = 0;
  (*minor_version) = 0;
  (*revision_version) = 0;
#endif
}

const char *
gribapiLibraryVersionString(void)
{
#ifdef HAVE_LIBGRIB_API
  if (!gribapi_libvers_init)
    {
      int major_version, minor_version, revision_version;
      gribapiLibraryVersion(&major_version, &minor_version, &revision_version);

      snprintf(gribapi_libvers, sizeof(gribapi_libvers), "%d.%d.%d", major_version, minor_version, revision_version);
      gribapi_libvers_init = true;
    }
#endif

  return gribapi_libvers;
}

void *
gribHandleNew(int editionNumber)
{
#ifdef HAVE_LIBGRIB_API
  grib_handle *gh = NULL;
  const char *fname = (editionNumber == 1) ? CDI_GRIB1_Template : CDI_GRIB2_Template;
  if (fname)
    {
      FILE *fp = fopen(fname, "r");
      if (fp)
        {
          int error;
          gh = grib_handle_new_from_file(NULL, fp, &error);
          fclose(fp);
          if (gh == NULL) Error("grib_handle_new_from_file failed!");
        }
      else
        {
          Error("Open failed on >%s<!", fname);
        }
    }

  if (gh == NULL)
    {
      gh = grib_handle_new_from_samples(NULL, (editionNumber == 1) ? "GRIB1" : "GRIB2");
      if (gh == NULL) Error("grib_handle_new_from_samples failed!");

      if (editionNumber == 1) GRIB_CHECK(my_grib_set_long(gh, "deleteLocalDefinition", 1L), 0);
      if (editionNumber == 2) GRIB_CHECK(my_grib_set_long(gh, "grib2LocalSectionPresent", 0L), 0);
      if (editionNumber == 2) GRIB_CHECK(my_grib_set_long(gh, "numberOfValues", 0L), 0);
    }

  return gh;
#else
  return NULL;
#endif
}

void
gribHandleDelete(void *gh)
{
#ifdef HAVE_LIBGRIB_API
  grib_handle_delete((struct grib_handle *) gh);
#endif
}

void
gribContainersNew(stream_t *streamptr)
{
  const int editionNumber = (streamptr->filetype == CDI_FILETYPE_GRB) ? 1 : 2;

#ifdef HAVE_LIBCGRIBEX
  if (editionNumber == 1 && !CDI_gribapi_grib1)
    {
    }
  else
#endif
    {
      const int nvars = streamptr->nvars;

#ifdef GRIBCONTAINER2D
      gribContainer_t **gribContainers;
      gribContainers = (gribContainer_t **) Malloc(nvars * sizeof(gribContainer_t *));

      for (int varID = 0; varID < nvars; ++varID)
        {
          const int nlevs = streamptr->vars[varID].nlevs;
          gribContainers[varID] = (gribContainer_t *) Malloc(nlevs * sizeof(gribContainer_t));

          for (int levelID = 0; levelID < nlevs; ++levelID)
            {
              gribContainers[varID][levelID].gribHandle = gribHandleNew(editionNumber);
              gribContainers[varID][levelID].init = false;
            }
        }

      streamptr->gribContainers = (void *) gribContainers;
#else
    gribContainer_t *gribContainers = (gribContainer_t *) Malloc((size_t) nvars * sizeof(gribContainer_t));

    for (int varID = 0; varID < nvars; ++varID)
      {
        gribContainers[varID].gribHandle = gribHandleNew(editionNumber);
        gribContainers[varID].init = false;
      }

    streamptr->gribContainers = (void *) gribContainers;
#endif
    }
}

void
gribContainersDelete(stream_t *streamptr)
{
  if (streamptr->gribContainers)
    {
      const int nvars = streamptr->nvars;

#ifdef GRIBCONTAINER2D
      gribContainer_t **gribContainers = (gribContainer_t **) streamptr->gribContainers;

      for (int varID = 0; varID < nvars; ++varID)
        {
          const int nlevs = streamptr->vars[varID].nlevs;
          for (int levelID = 0; levelID < nlevs; ++levelID)
            {
              gribHandleDelete(gribContainers[varID][levelID].gribHandle);
            }
          Free(gribContainers[varID]);
        }
#else
      gribContainer_t *gribContainers = (gribContainer_t *) streamptr->gribContainers;

      for (int varID = 0; varID < nvars; ++varID)
        {
          gribHandleDelete(gribContainers[varID].gribHandle);
        }
#endif

      Free(gribContainers);

      streamptr->gribContainers = NULL;
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
