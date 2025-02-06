#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdi.h"

// A version string

#ifdef VERSION
static const char cdi_libvers[] = VERSION;
#else
#error "VERSION undefined"
#endif

const char *
cdiLibraryVersion(void)
{
  return cdi_libvers;
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
